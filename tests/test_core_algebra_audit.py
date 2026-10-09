"""Independent Phase 4G-1 references; production algorithms are unchanged.

Random datasets use explicit seeds. Fraction elimination, Decimal norms,
Leibniz determinants, Rodrigues rotations and analytic camera equations are
oracles, rather than round trips alone. Confirmed defects stay strict xfails.
"""
import copy
import ctypes
from decimal import Decimal, localcontext
from fractions import Fraction
import itertools
import math
import random

import pytest

from gem import common, matrix, quaternion, vector
from .helpers import inverse as rational_inverse


def V(values):
    return vector.Vector(len(values), list(values))


def decimal_direction(values):
    with localcontext() as context:
        context.prec = 800
        entries = [Decimal.from_float(float(x)) for x in values]
        length = sum(x*x for x in entries).sqrt()
        direction = [float(x/length) for x in entries] if length else [0.0]*len(entries)
        return float(length), direction


def leibniz(rows):
    """Exact permutation expansion, independent of gem's cofactor plumbing."""
    total = Fraction(0)
    for permutation in itertools.permutations(range(len(rows))):
        inversions = sum(permutation[i] > permutation[j]
                         for i in range(len(rows)) for j in range(i+1, len(rows)))
        term = Fraction((-1)**inversions)
        for i, j in enumerate(permutation):
            term *= Fraction(rows[i][j])
        total += term
    return total


def rational_product(a, b):
    return [[float(sum(Fraction(a[i][k])*Fraction(b[k][j])
                       for k in range(len(a))))
             for j in range(len(a))] for i in range(len(a))]


def hamilton(a, b):
    """Scalar/vector decomposition with exact rational operations."""
    aw, ax, ay, az = map(Fraction, a)
    bw, bx, by, bz = map(Fraction, b)
    return [aw*bw - ax*bx - ay*by - az*bz,
            aw*bx + bw*ax + ay*bz - az*by,
            aw*by + bw*ay + az*bx - ax*bz,
            aw*bz + bw*az + ax*by - ay*bx]


def rodrigues(axis, angle, point):
    """Active column-vector Rodrigues formula, not a quaternion sandwich."""
    length = math.sqrt(sum(x*x for x in axis))
    nx, ny, nz = [x/length for x in axis]
    x, y, z = point
    cross = [ny*z-nz*y, nz*x-nx*z, nx*y-ny*x]
    parallel = nx*x + ny*y + nz*z
    c, s = math.cos(angle), math.sin(angle)
    return [p*c + perpendicular*s + n*parallel*(1-c)
            for p, perpendicular, n in zip(point, cross, [nx, ny, nz])]


def assert_rows(actual, expected, rel=3e-14, abs=3e-14):
    assert len(actual) == len(expected)
    for row, reference in zip(actual, expected):
        assert row == pytest.approx(reference, rel=rel, abs=abs)


@pytest.mark.defect('4G1-A01: subnormal imaginary-axis scaling in quaternion powers')
@pytest.mark.parametrize('components', [(1, 1, 0), (1, -1, 1)])
def test_subnormal_quaternion_square_root_axis(components):
    tiny = math.ldexp(1.0, -1074)
    data = [-1.0] + [component*tiny for component in components]
    _, axis = decimal_direction(data[1:])
    original = data[:]
    q = quaternion.Quaternion(data)
    # Existing unit-input tests accept +/-1 with tiny imaginary components.
    # Their exact norm differs from 1 by less than a binary64 rounding unit.
    result = q.pow(0.5)
    assert result is not q and result.data is not data
    assert data == original and q.data is data
    assert result.data == pytest.approx([0.0] + axis, rel=2e-15, abs=1e-15)
    assert math.hypot(*result.data) == pytest.approx(1.0, abs=2e-15)


@pytest.mark.defect('4G1-A01: subnormal imaginary-axis scaling in quaternion logarithms')
@pytest.mark.parametrize('components', [(1, 1, 0), (1, -1, 1)])
def test_subnormal_quaternion_logarithm_axis(components):
    tiny = math.ldexp(1.0, -1074)
    data = [-1.0] + [component*tiny for component in components]
    _, axis = decimal_direction(data[1:])
    original = data[:]
    q = quaternion.Quaternion(data)
    result = q.log()
    assert type(result) is list and result is not data
    assert data == original and q.data is data
    assert result == pytest.approx([0.0] + [math.pi*x for x in axis], rel=2e-15, abs=0)


@pytest.mark.defect('4G1-A02: large integer powers lose the algebraic rotation phase')
@pytest.mark.parametrize('data', [[0.0, 1.0, 0.0, 0.0], [0.5, 0.5, 0.5, 0.5]])
@pytest.mark.parametrize('exponent', [10**16, -10**16])
def test_large_integer_quaternion_power(data, exponent):
    # Independently repeated exact Hamilton products: i^4=1; this .5 control
    # has order six. Both data sets have exactly unit norm, not just near-unit.
    period = 4 if data[0] == 0 else 6
    expected = [Fraction(1), Fraction(0), Fraction(0), Fraction(0)]
    for _ in range(exponent % period):
        expected = hamilton(expected, data)
    original = data[:]
    q = quaternion.Quaternion(data)
    result = q.pow(exponent)
    assert q.data is data and data == original
    assert result is not q and result.data is not data
    # Check represented rotation too: this is not merely a sign mismatch.
    assert abs(sum(float(x)*y for x, y in zip(expected, result.data))) == pytest.approx(1, abs=2e-14)
    assert result.data == pytest.approx([float(x) for x in expected], abs=2e-14)


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('seed', range(8))
def test_vector_exact_dyadic_algebra_and_aliasing(size, seed):
    rng = random.Random(470100 + seed)
    a = [rng.randint(-40, 40)/8.0 for _ in range(size)]
    b = [rng.randint(-40, 40)/8.0 for _ in range(size)]
    va, vb = vector.Vector(size, a), vector.Vector(size, b)
    original = (a[:], b[:])
    assert va.dot(vb) == float(sum(Fraction(x)*Fraction(y) for x, y in zip(a, b)))
    for result, expected in [(va+vb, [x+y for x, y in zip(a, b)]),
                             (va-vb, [x-y for x, y in zip(a, b)]),
                             (va*2, [x*2 for x in a])]:
        assert result.vector == expected and result.vector is not a
        result.vector[0] = 999
        assert (a, b) == original
    alias = vector.Vector(size, a)
    assert va.__iadd__(va) is va
    assert va.vector == [x*2 for x in original[0]]
    assert alias.vector is a and a == original[0]
    assert (alias == V(a)) and not (alias != V(a))


@pytest.mark.parametrize('values', [
    [math.ldexp(1.0, -1074)]*2,
    [1e-300, -2e-300, 3e-300],
    [1e300, -1e-300, 2e300],
    [1e308, -1e308, 1e308, -1e308],
    [0.0, -0.0, 0.0],
])
def test_stable_vector_norm_decimal_and_viewport(values):
    expected_length, expected_direction = decimal_direction(values)
    v = V(values)
    storage = v.vector
    normalized = v.normalize()
    assert v.magnitude() == pytest.approx(expected_length, rel=3e-15, abs=0)
    assert normalized.vector == pytest.approx(expected_direction, rel=3e-15, abs=math.ldexp(1.0, -1074))
    assert normalized.vector is not storage and storage == values
    if any(values):
        expected = [(expected_direction[0]+1)*160 + values[0],
                    (expected_direction[1]+1)*90 + values[1], 320, 180]
        assert common.getViewPort(v, 320, 180) == pytest.approx(expected, rel=3e-15, abs=1e-13)
    else:
        with pytest.raises(ZeroDivisionError):
            common.getViewPort(v, 320, 180)
    assert v.vector is storage and storage == values


def test_cross_reflection_and_barycentric_analytic_answers():
    a, b = V([2, -3, 5]), V([-7, 11, 13])
    cross = vector.cross(a, b)
    assert cross.vector == [-94, -61, 1]
    assert cross.dot(a) == cross.dot(b) == 0
    normal = V([0, 1, 0])
    assert vector.reflect(a, normal).vector == [2, 3, 5]
    # Non-orthogonal triangle: p=.25*a+.5*b+.25*c, independently constructed.
    va, vb, vc = V([1, 2, 3]), V([5, 2, 3]), V([3, 6, 3])
    assert V([3.5, 3, 3]).barycentric(va, vb, vc) == [0.25, 0.5, 0.25]
    assert a.vector == [2, -3, 5] and b.vector == [-7, 11, 13]


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('seed', range(8))
def test_matrix_algebra_independent_exact_references(size, seed):
    rng = random.Random(470200 + seed)
    rows = [[float(rng.randint(-3, 3) + (12 if i == j else 0))
             for j in range(size)] for i in range(size)]
    others = [[[float(rng.randint(-3, 3)) for _ in range(size)] for _ in range(size)] for _ in range(2)]
    original = copy.deepcopy(rows)
    a, b, c = [matrix.Matrix(size, x) for x in [rows] + others]
    assert a.det() == float(leibniz(rows))
    expected_inverse = rational_inverse(rows)
    assert_rows(a.inverse().matrix, expected_inverse)
    assert (a*b).matrix == rational_product(rows, others[0])
    assert ((a*b)*c).matrix == (a*(b*c)).matrix  # Exact small integer products.
    for product in [rational_product(rows, a.inverse().matrix),
                    rational_product(a.inverse().matrix, rows)]:
        assert_rows(product, matrix.identity(size))
    point = [float(rng.randint(-3, 3)) for _ in range(size)]
    expected = [sum(point[i]*rows[i][j] for i in range(size)) for j in range(size)]
    assert (a*V(point)).vector == expected
    assert vector.transform(size, point, rows) == expected
    assert a.__imul__(a) is a
    assert a.matrix == rational_product(original, original)
    assert rows == original
    assert [list(row) for row in a.c_matrix] == [[ctypes.c_float(x).value for x in row] for row in a.matrix]


@pytest.mark.parametrize('size', [3, 4])
@pytest.mark.parametrize('exponent', [-996, -500, 0, 500, 996])
def test_nonsymmetric_scaled_inverse_fraction_reference(size, exponent):
    # Fixed ordinary condition number; only a uniform power-of-two scale varies.
    base = [[4 if i == j else -1 if j == i+1 else 1 if i == j+1 else 0
             for j in range(size)] for i in range(size)]
    rows = [[math.ldexp(float(x), exponent) for x in row] for row in base]
    original = copy.deepcopy(rows)
    result = getattr(matrix, 'inverse'+str(size))(rows)
    assert_rows(result, rational_inverse(rows), abs=0)
    for product in [rational_product(rows, result), rational_product(result, rows)]:
        assert_rows(product, matrix.identity(size))
    assert rows == original
    rows[-1] = rows[0][:]
    assert leibniz(rows) == 0
    with pytest.raises(ZeroDivisionError):
        getattr(matrix, 'inverse'+str(size))(rows)


@pytest.mark.parametrize('seed', range(12))
def test_quaternion_matrix_rotation_against_rodrigues(seed):
    rng = random.Random(470300 + seed)
    axis = [rng.uniform(-2, 2) for _ in range(3)]
    angle = rng.uniform(-math.pi, math.pi)
    point = [rng.uniform(-4, 4) for _ in range(3)]
    length = math.sqrt(sum(x*x for x in axis))
    data = [math.cos(angle/2)] + [x/length*math.sin(angle/2) for x in axis]
    q = quaternion.Quaternion(data)
    before = data[:]
    expected = rodrigues(axis, angle, point)
    assert quaternion.quat_from_axis_angle(axis, math.degrees(angle)).data == pytest.approx(data, abs=2e-15)
    rotated = quaternion.quat_rotate_vector(q, V(point))
    assert rotated.vector == pytest.approx(expected, abs=4e-15)
    assert (q.toMatrix()*V(point+[0])).vector == pytest.approx(expected+[0], abs=4e-15)
    assert (matrix.Matrix(3, matrix.rotate3(axis, math.degrees(angle)))*V(point)).vector == pytest.approx(expected, abs=4e-15)
    assert math.hypot(*rotated.vector) == pytest.approx(math.hypot(*point), rel=3e-15)
    back = quaternion.quat_from_matrix(q.toMatrix())
    assert abs(sum(a*b for a, b in zip(back.data, data))) == pytest.approx(1, abs=2e-15)
    assert quaternion.quat_rotate_vector(q.negate(), V(point)).vector == rotated.vector
    assert q.data is data and data == before


def test_quaternion_noncommuting_composition_and_inverse():
    # Exact unit quaternions for 120-degree turns about different diagonal axes.
    a = [0.5, 0.5, 0.5, 0.5]
    b = [0.5, -0.5, 0.5, -0.5]
    qa, qb = quaternion.Quaternion(a), quaternion.Quaternion(b)
    expected = [float(x) for x in hamilton(a, b)]
    assert (qa*qb).data == expected and (qb*qa).data != expected
    assert (qa*qb).toMatrix().matrix == (qb.toMatrix()*qa.toMatrix()).matrix
    for data in [[2, -3, 4, -5], a, b]:
        q = quaternion.Quaternion(list(data))
        norm_squared = sum(Fraction(x)**2 for x in data)
        inv = [float(Fraction(x)*(1 if i == 0 else -1)/norm_squared) for i, x in enumerate(data)]
        assert q.inverse().data == inv
        for product in [q*q.inverse(), q.inverse()*q]:
            assert product.data == pytest.approx([1, 0, 0, 0], abs=2e-15)


@pytest.mark.parametrize('t', [0.0, 0.125, 0.5, 0.875, 1.0])
@pytest.mark.parametrize('sign', [1, -1])
def test_slerp_and_squad4_independent_axis_angles(t, sign):
    def Q(angle):
        return quaternion.Quaternion([math.cos(angle), 0, 0, math.sin(angle)])
    # Arguments are quaternion half-angles, not rotation angles.
    a, b, s0, s1 = Q(-0.3), Q(0.6), Q(0.1), Q(0.4)
    data = [q.data[:] for q in [a, b, s0, s1]]
    end = b if sign == 1 else b.negate()
    expected_angle = -0.3*(1-t) + 0.6*t
    assert a.slerp(end, t).data == pytest.approx(Q(expected_angle).data, abs=2e-15)
    blend = 2*t*(1-t)
    squad_angle = (1-blend)*expected_angle + blend*(0.1*(1-t)+0.4*t)
    output = quaternion.squad4(a, end, s0, s1, t)
    assert output.data == pytest.approx(Q(squad_angle).data, abs=2e-15)
    assert math.hypot(*output.data) == pytest.approx(1, abs=2e-15)
    assert [q.data for q in [a, b, s0, s1]] == data
    assert all(output.data is not q.data for q in [a, b, s0, s1])


def test_legacy_squad_midpoint_on_long_arc():
    controls = [quaternion.Quaternion([math.cos(angle), 0, 0, math.sin(angle)])
                for angle in [0, math.pi/3, 2*math.pi/3]]
    # At t=.5 all three no-invert blends use their spherical branches:
    # a half-angle pi/3, b pi/6, final pi/4 (90-degree orientation).
    original = [q.data[:] for q in controls]
    output = controls[0].squad(controls[1], controls[2], 0.5)
    assert output.data == pytest.approx([math.sqrt(0.5), 0, 0, math.sqrt(0.5)], abs=2e-15)
    assert [q.data for q in controls] == original


@pytest.mark.parametrize('axis_index', [0, 1, 2])
@pytest.mark.parametrize('sign', [1, -1])
def test_exact_half_turn_conversion_and_slerp_tie(axis_index, sign):
    data = [0.0, 0.0, 0.0, 0.0]
    data[axis_index+1] = float(sign)
    q = quaternion.Quaternion(data)
    expected_rows = matrix.identity(4)
    for i in range(3):
        expected_rows[i][i] = 1 if i == axis_index else -1
    assert q.toMatrix().matrix == expected_rows
    back = quaternion.quat_from_matrix(matrix.Matrix(4, expected_rows))
    assert abs(sum(x*y for x, y in zip(back.data, data))) == 1
    expected = [math.sqrt(0.5), 0.0, 0.0, 0.0]
    expected[axis_index+1] = sign*math.sqrt(0.5)
    # Exactly dot==0: the documented tie retains the supplied endpoint sign.
    assert quaternion.Quaternion().slerp(q, 0.5).data == pytest.approx(expected, abs=2e-15)


@pytest.mark.parametrize('size', [2, 3, 4])
def test_singular_inplace_inverse_preserves_receiver_and_export(size):
    rows = matrix.identity(size)
    rows[-1] = rows[0][:]
    m = matrix.Matrix(size, rows)
    before = (copy.deepcopy(rows), [list(row) for row in m.c_matrix])
    for method in [m.inverse, m.i_inverse]:
        with pytest.raises(ZeroDivisionError):
            method()
        assert m.matrix is rows
        assert (rows, [list(row) for row in m.c_matrix]) == before


def test_projection_behind_camera_preserves_unclamped_depth():
    p = matrix.perspective(90, 1, 1, 9)
    # For [1,2,+2,1], NDC=[-.5,-1,2.375]; negative W is not clipped.
    result = matrix.project(V([1, 2, 2, 1]), matrix.Matrix(4), p, [0, 0, 100, 100])
    assert result.vector == pytest.approx([25, 0, 27/16], abs=2e-14)
    assert matrix.unproject(25, 0, 27/16, matrix.Matrix(4), p, [0, 0, 100, 100]).vector == pytest.approx([1, 2, 2], abs=2e-14)


@pytest.mark.parametrize('kind', ['perspective', 'orthographic'])
@pytest.mark.parametrize('seed', range(8))
def test_camera_noncommuting_independent_window_reference(kind, seed):
    rng = random.Random(470400 + seed)
    x, y, z = rng.uniform(-1, 1), rng.uniform(-1, 1), rng.uniform(-8, -3)
    # Analytic +Z quarter-turn, then translation; this does not commute with P.
    rows = [[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 1, 0], [2, -1, -2, 1]]
    ex, ey, ez = 2-y, x-1, z-2
    viewport = [-13, 21, 420, 170]
    if kind == 'perspective':
        projection = matrix.perspective(60, 2.5, 1, 20)
        f = 1/math.tan(math.pi/6)
        nx, ny = f*ex/(-2.5*ez), f*ey/(-ez)
        depth = 20/19 + 20/(19*ez)
    else:
        projection = matrix.orthographic(-4, 6, -3, 5, 1, 20)
        nx, ny = 2*(ex+4)/10-1, 2*(ey+3)/8-1
        depth = (-ez-1)/19
    window = [-13+(nx+1)*210, 21+(ny+1)*85, depth]
    model = matrix.Matrix(4, rows) if seed % 2 else rows
    proj = projection if seed % 3 else projection.matrix
    original = copy.deepcopy((rows, projection.matrix, viewport))
    assert matrix.project(V([x, y, z, 1]), model, proj, viewport).vector == pytest.approx(window, abs=1e-12)
    assert matrix.unproject(*window, model, proj, viewport).vector == pytest.approx([x, y, z], abs=2e-13)
    assert (rows, projection.matrix, viewport) == original


def test_degenerate_and_unsupported_contracts():
    ident = matrix.Matrix(4)
    with pytest.raises(ZeroDivisionError):
        matrix.project(V([1, 2, 3, 0]), ident, ident, [0, 0, 100, 100])
    with pytest.raises(IndexError):
        ident*V([1, 2, 3])  # No general implicit homogeneous promotion.
    with pytest.raises(ZeroDivisionError):
        matrix.lookAt(V([1, 2, 3]), V([1, 2, 3]), V([0, 1, 0]))
    with pytest.raises(ZeroDivisionError):
        matrix.lookAt(V([0, 0, 0]), V([0, 0, -1]), V([0, 0, 1]))
    with pytest.raises(ValueError):
        quaternion.Quaternion([-1, 0, 0, 0]).log()
    assert V([1, 2]).__add__(object()) is NotImplemented
    assert V([1, 2]).__eq__(object()) is NotImplemented
    assert ident.__mul__(object()) is NotImplemented
    assert quaternion.Quaternion().__mul__(2) is NotImplemented


def test_angle_helpers_independent_periods_and_keyword_compatibility():
    assert common.radiansToDegrees(degrees=math.pi/3) == pytest.approx(60)
    assert common.degreesToRadians(radians=-270) == pytest.approx(-3*math.pi/2)
    assert vector.toAngle([-1, 0]) == math.pi
    assert vector.lperp([3, -5]).vector == [5, 3]
    assert vector.rperp([3, -5]).vector == [-5, -3]
