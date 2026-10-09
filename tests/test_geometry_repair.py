"""Geometric identities independent of refraction and Newell implementations."""
from decimal import Decimal, localcontext
from fractions import Fraction
import math

import pytest

from gem import plane, vector


def V(values):
    return vector.Vector(len(values), list(values))


def area_reference(points):
    """Exact triangle-fan area, normalized using high-precision arithmetic."""
    anchor = [Fraction(x) for x in points[0]]
    area = [Fraction(0)] * 3
    for i in range(1, len(points) - 1):
        a = [Fraction(x) - y for x, y in zip(points[i], anchor)]
        b = [Fraction(x) - y for x, y in zip(points[i + 1], anchor)]
        cross = [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
                 a[0]*b[1] - a[1]*b[0]]
        area = [x + y for x, y in zip(area, cross)]
    with localcontext() as context:
        context.prec = 110
        values = [Decimal(x.numerator) / Decimal(x.denominator) for x in area]
        length = sum(x*x for x in values).sqrt()
        return [float(x / length) for x in values]


@pytest.mark.parametrize('size', [2, 3, 4])
@pytest.mark.parametrize('side', [-1, 1])
@pytest.mark.parametrize('tangent_sign', [-1, 1])
@pytest.mark.parametrize('sine', [0., 5e-324, 1e-310, 1e-100, 1e-9, .6, 1.])
def test_equal_indices_preserve_direction_at_both_interfaces(size, side, tangent_sign, sine):
    incident = V([tangent_sign * math.sqrt(1. - sine*sine), -side*sine]
                 + [0.] * (size - 2))
    normal = V([0., float(side)] + [0.] * (size - 2))
    incident_storage, normal_storage = incident.vector, normal.vector
    expected = incident.vector[:]
    # Snell's law at n1=n2 leaves the direction unchanged, including subnormals.
    result = vector.refract(1., incident, normal)
    assert result.vector == expected
    assert result is not incident and result is not normal
    assert result.vector is not incident_storage and result.vector is not normal_storage
    assert incident.vector is incident_storage and incident.vector == expected
    assert normal.vector is normal_storage


@pytest.mark.parametrize('side', [-1, 1])
def test_equal_indices_non_axis_aligned_interface(side):
    normal = V([.6*side, .8*side, 0.])
    incident = V([.8*side, -.6*side, 0.])
    result = vector.refract(1., incident, normal)
    assert result.vector == incident.vector


def test_equal_index_correction_keeps_existing_operand_dispatch():
    with pytest.raises(TypeError):
        vector.refract(Fraction(1), V([1., -1e-9]), V([0., 1.]))
    # No automatic normal flipping: this is outside the opposing-normal contract.
    assert vector.refract(1., V([.8, .6]), V([0., 1.])).vector == pytest.approx([.8, -.6])


@pytest.mark.parametrize('eta', [.5, 2./3., 1.5, 2.])
@pytest.mark.parametrize('side', [-1, 1])
def test_ordinary_snell_and_total_internal_reflection(eta, side):
    # Tangential sine scales by eta; transmitted normal component is negative.
    for sine in [.2, .6, .9]:
        incident = V([sine, -side*math.sqrt(1. - sine*sine), 0.])
        result = vector.refract(eta, incident, V([0., float(side), 0.]))
        transmitted = eta*sine
        expected = ([0., 0., 0.] if transmitted > 1. else
                    [transmitted, -side*math.sqrt(1. - transmitted*transmitted), 0.])
        assert result.vector == pytest.approx(expected, rel=3e-15, abs=3e-15)


POLYGONS = [
    [[0, 0, 0], [2, 0, 4], [0, 2, -6]],
    [[0, 0, 7], [4, 0, 15], [4, 4, 3], [0, 4, -5]],
    # Concave planar polygon, z=2x-3y+7.
    [[0, 0, 7], [4, 0, 15], [4, 4, 3], [2, 2, 5], [0, 4, -5]],
    # Nonplanar Newell area remains an area-weighted sum, not a plane fit.
    [[0, 0, 0], [4, 0, 2], [4, 4, 5], [0, 4, -1]],
]


@pytest.mark.parametrize('local', POLYGONS)
@pytest.mark.parametrize('offset', [(0., 0., 0.), (2.**52, -2.**52, 2.**51),
                                   (-2.**52, 2.**51, -2.**52)])
@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('closed', [False, True])
def test_polygon_area_is_translation_and_start_vertex_invariant(local, offset, reverse, closed):
    points = [[x+y for x, y in zip(point, offset)] for point in local]
    if reverse:
        points.reverse()
    expected = area_reference(points)
    for start in range(len(points)):
        ordered = points[start:] + points[:start]
        if closed:
            ordered = ordered + [ordered[0][:]]
        vertices = [V(point) for point in ordered]
        storages = [point.vector for point in vertices]
        result = plane.Plane().bestFitNormal(vertices)
        assert result.vector == pytest.approx(expected, rel=3e-15, abs=3e-15)
        assert math.hypot(*result.vector) == pytest.approx(1., abs=3e-15)
        for vertex, storage, original in zip(vertices, storages, ordered):
            assert vertex.vector is storage and vertex.vector == original
            assert result.vector is not storage


@pytest.mark.parametrize('exponent', [-200, 0, 200])
def test_polygon_uniform_scale_retains_orientation(exponent):
    scale = math.ldexp(1., exponent)
    points = [[x*scale for x in point] for point in POLYGONS[0]]
    assert plane.Plane().bestFitNormal([V(point) for point in points]).vector == pytest.approx(
        area_reference(points), rel=3e-15, abs=3e-15)


@pytest.mark.parametrize('offset,edge', [(1e100, 1e90), (1e150, 1e140),
                                        (1e-100, 1e-110), (1e-140, 1e-150)])
@pytest.mark.parametrize('reverse', [False, True])
def test_polygon_offsets_across_finite_exponents(offset, edge, reverse):
    points = [[offset+x*edge for x in point] for point in POLYGONS[1]]
    if reverse:
        points.reverse()
    # Reference uses the exact represented coordinates, not ideal pre-rounding data.
    assert plane.Plane().bestFitNormal([V(point) for point in points]).vector == pytest.approx(
        area_reference(points), rel=3e-15, abs=3e-15)


@pytest.mark.parametrize('points', [[], [[0, 0, 0]], [[1, 2, 3], [4, 5, 6]],
                                   [[2.**52, 0., 0.]] * 3,
                                   [[0, 0, 0], [1, 2, 3], [2, 4, 6]],
                                   [[0, 0, 0], [2, 2, 0], [0, 2, 0], [2, 0, 0]]])
def test_degenerate_polygon_keeps_zero_division_error(points):
    with pytest.raises(ZeroDivisionError):
        plane.Plane().bestFitNormal([V(point) for point in points])
