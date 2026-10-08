"""Independent geometry regressions for Phase 2D-1."""
import math
import random
import pytest
from gem import plane, vector


def V(values):
    return vector.Vector(len(values), list(values))


def coefficients(p):
    return [p.a, p.b, p.c, p.d]


@pytest.mark.parametrize('points,expected', [
    ([(0,2,0),(0,2,1),(1,2,0)], [0,1,0,-2]),
    ([(5,0,0),(5,1,0),(5,0,1)], [1,0,0,-5]),
    ([(0,0,-3),(1,0,-3),(0,1,-3)], [0,0,1,3]),
    ([(1,2,3),(4,0,3),(13,20,-10)], [2/7,3/7,6/7,-26/7]),
])
def test_three_point_known_answers(points, expected):
    inputs = [V(point) for point in points]
    p = plane.Plane()
    assert p.fromPoints(*inputs) is None
    assert coefficients(p) == pytest.approx(expected, abs=1e-14)
    assert p.normal.vector == pytest.approx(expected[:3], abs=1e-14)
    for point in points:
        assert p.dot(V(list(point)+[1])) == pytest.approx(0, abs=1e-13)
    for value,original in zip(inputs,points):
        assert value.vector == list(original)
    reversed_plane = plane.Plane()
    reversed_plane.fromPoints(inputs[0], inputs[2], inputs[1])
    assert coefficients(reversed_plane) == pytest.approx([-x for x in expected], abs=1e-14)
    inputs[0].vector[0] += 99
    assert coefficients(p) == pytest.approx(expected, abs=1e-14)


@pytest.mark.parametrize('scale', [-4.0, -0.5, 0.125, 1.0, 3.0])
def test_normalization_scale_invariance(scale):
    # Normal length is exactly 7; every point satisfies 2*x+3*y+6*z=26.
    values = [2*scale,3*scale,6*scale,-26*scale]
    expected = [x/(7*abs(scale)) for x in values]
    assert plane.normalize(values) == pytest.approx(expected, abs=1e-14)
    assert values == [2*scale,3*scale,6*scale,-26*scale]
    p = plane.Plane()
    p.fromCoeffs(*values)
    q = p.normalize()
    assert q is not p
    assert coefficients(p) == values
    assert coefficients(q) == pytest.approx(expected, abs=1e-14)
    assert q.normal.vector == pytest.approx(expected[:3], abs=1e-14)
    assert q.normal.magnitude() == pytest.approx(1, abs=1e-14)
    for point in [(1,2,3),(4,0,3),(13,20,-10)]:
        assert p.dot(V(list(point)+[1])) == pytest.approx(0, abs=1e-12)
        assert q.dot(V(list(point)+[1])) == pytest.approx(0, abs=1e-13)
    assert p.i_normalize() is p
    assert coefficients(p) == pytest.approx(expected, abs=1e-14)
    assert p.normal.vector == pytest.approx(expected[:3], abs=1e-14)
    assert coefficients(p.normalize()) == pytest.approx(expected, abs=1e-14)


@pytest.mark.parametrize('seed', range(20))
def test_analytic_graph_plane_incidence_and_distance(seed):
    rng = random.Random(seed)
    a,b,c = [rng.uniform(-5,5) for _ in range(3)]
    # z=a*x+b*y+c gives an analytic normal independent of gem cross products.
    length = math.sqrt(a*a+b*b+1)
    expected = [-a/length,-b/length,1/length,-c/length]
    p = plane.Plane()
    p.fromPoints(V([0,0,c]), V([1,0,a+c]), V([0,1,b+c]))
    assert coefficients(p) == pytest.approx(expected, abs=1e-13)
    for _ in range(3):
        x,y = rng.uniform(-10,10), rng.uniform(-10,10)
        z = a*x+b*y+c
        assert p.dot(V([x,y,z,1])) == pytest.approx(0, abs=1e-12)
        offset = rng.uniform(-5,5)
        assert p.dot(V([x,y,z+offset,1])) == pytest.approx(offset/length, abs=1e-12)
    q = plane.Plane()
    q.fromCoeffs(-a*3,-b*3,3,-c*3)
    assert coefficients(q.normalize()) == pytest.approx(expected, abs=1e-13)


@pytest.mark.parametrize('values,point', [
    ([0,2,0,-4], [0,2,0]),
    ([0,0,3,9], [0,0,-3]),
    ([2,3,6,-26], [1,2,3]),
    ([-2,-3,-6,26], [1,2,3]),
])
def test_coefficient_scale_normal_and_side(values, point):
    p = plane.Plane()
    assert p.fromCoeffs(*values) is None
    assert coefficients(p) == values
    assert p.normal.vector == values[:3]
    assert p.dot(V(point+[1])) == 0
    for sign in [-1,1]:
        off_plane = [x+sign*n for x,n in zip(point,values[:3])]
        assert p.point_location(p, off_plane) == sign
        assert p.dot(V(off_plane+[1])) == sign*sum(n*n for n in values[:3])
    for q in [p.clone(), p.flip()]:
        assert q.normal is not p.normal
        assert q.normal.vector == [q.a,q.b,q.c]
        assert q.dot(V(point+[1])) == 0
    before = coefficients(p)
    assert p.i_flip() is p
    assert coefficients(p) == [-x for x in before]
    assert p.normal.vector == [p.a,p.b,p.c]


@pytest.mark.parametrize('points,normal,offset', [
    ([(0,0,2),(2,0,2),(2,1,2),(0,1,2)], [0,0,1], 2),
    ([(0,0,-3),(2,0,-3),(2,1,-3),(0,1,-3)], [0,0,1], -3),
    ([(0,2,0),(0,2,1),(1,2,1),(1,2,0)], [0,1,0], 2),
    ([(1,2,3),(4,0,3),(16,18,-10),(13,20,-10)], [2/7,3/7,6/7], 26/7),
])
@pytest.mark.parametrize('closed', [False, True], ids=['implicit-wrap','repeated-first'])
def test_polygon_plane_wrapping_incidence_and_winding(points, normal, offset, closed):
    vertices = [V(point) for point in points]
    if closed:
        vertices.append(vertices[0].clone())
    originals = [list(v.vector) for v in vertices]
    helper = plane.Plane()
    for shift in range(len(points)):
        rotated = vertices[shift:len(points)]+vertices[:shift]
        if closed:
            rotated.append(rotated[0].clone())
        n = helper.bestFitNormal(rotated)
        D = helper.bestFitD(rotated,n)
        assert n.vector == pytest.approx(normal, abs=1e-14)
        assert D == pytest.approx(offset, abs=1e-14)
        p = plane.Plane()
        p.fromCoeffs(n.vector[0], n.vector[1], n.vector[2], -D)
        for point in points:
            assert p.dot(V(list(point)+[1])) == pytest.approx(0, abs=1e-13)
    reversed_vertices = list(reversed(vertices))
    n = helper.bestFitNormal(reversed_vertices)
    assert n.vector == pytest.approx([-x for x in normal], abs=1e-14)
    assert helper.bestFitD(reversed_vertices,n) == pytest.approx(-offset, abs=1e-14)
    assert [v.vector for v in vertices] == originals
    assert coefficients(helper) == [0,0,0,0]


def test_best_fit_offset_retains_signed_mean_and_normal_scale():
    p = plane.Plane()
    vertices = [V([0,0,-2]), V([1,0,-4]), V([0,1,-6])]
    assert p.bestFitD(vertices,V([0,0,1])) == -4
    assert p.bestFitD(vertices,V([0,0,2])) == -8
    assert p.bestFitD(vertices+[vertices[0]],V([0,0,1])) == -3.5


def test_normalization_repairs_stale_normal_from_scalar_fields():
    p = plane.Plane()
    p.a,p.b,p.c,p.d = 0,0,3,9
    p.normal = V([99,0,0])
    q = p.normalize()
    assert coefficients(q) == [0,0,1,3]
    assert q.normal.vector == [0,0,1]
    assert p.normal.vector == [99,0,0]
    assert p.i_normalize() is p
    assert p.normal.vector == [0,0,1]


def test_existing_degenerate_errors():
    # Characterize current failures without adopting a broader domain policy.
    with pytest.raises(ZeroDivisionError):
        plane.normalize([0,0,0,1])
    with pytest.raises(ZeroDivisionError):
        plane.Plane().fromPoints(V([0,0,0]),V([1,1,1]),V([2,2,2]))
    with pytest.raises(ZeroDivisionError):
        plane.Plane().bestFitNormal([])
    with pytest.raises(ZeroDivisionError):
        plane.Plane().bestFitD([],V([0,0,1]))
