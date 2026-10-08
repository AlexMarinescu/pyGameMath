import pytest
from gem import bezier, vector


def V(*xs):
    return vector.Vector(len(xs), list(xs))


@pytest.mark.parametrize('t',[0,0.25,0.5,0.75,1])
def test_quadratic_bezier_scalar(t):
    assert bezier.quadraticBezierPoint(t,0,1,2) == pytest.approx(2*t)


@pytest.mark.parametrize('t',[0,0.25,0.5,0.75,1])
def test_cubic_bezier_scalar(t):
    assert bezier.cubicBezierPoint(t,0,1,2,3) == pytest.approx(3*t)


def test_bezier_vectors():
    assert bezier.quadraticBezierPoint(0.5,V(0,0),V(1,1),V(2,0)).vector == [1,0.5]


def test_bezier_path_count_python3():
    path = bezier.BezierPath()
    path.setControlPoints([0,1,2,3])
    assert isinstance(path.curveCount,int)
    assert path.curveCount == 1
    assert path.calculateBezerPoint(0, 0.5) == pytest.approx(1.5)




def casteljau(t, coordinates):
    points = list(coordinates)
    while len(points) > 1:
        points = [(1-t)*a+t*b for a,b in zip(points, points[1:])]
    return points[0]


@pytest.mark.parametrize('degree', [2, 3])
@pytest.mark.parametrize('dimension', [2, 3, 4])
@pytest.mark.parametrize('t', [-0.5, 0, 0.125, 0.5, 0.875, 1, 1.5])
def test_vector_evaluation_independent_reference(degree, dimension, t):
    coordinates = [[(-1)**i * (i+1) * (j+2) + i*i-j for j in range(dimension)]
                   for i in range(degree+1)]
    controls = [V(*p) for p in coordinates]
    storage = [p.vector for p in controls]
    evaluate = bezier.quadraticBezierPoint if degree == 2 else bezier.cubicBezierPoint
    result = evaluate(t, *controls)
    expected = [casteljau(t, [p[j] for p in coordinates]) for j in range(dimension)]
    assert result.vector == pytest.approx(expected)
    assert [p.vector for p in controls] == coordinates
    assert all(p.vector is saved for p,saved in zip(controls, storage))
    assert all(result is not p and result.vector is not p.vector for p in controls)
    assert evaluate(1-t, *reversed(controls)).vector == pytest.approx(expected)


@pytest.mark.parametrize('t', [0, 0.2, 0.5, 0.9, 1])
def test_degree_elevation(t):
    # Quadratic controls [-2, 4, 1] elevate to [-2, 2, 3, 1].
    expected = -2*(1-t)**2 + 8*(1-t)*t + t*t
    assert bezier.quadraticBezierPoint(t, -2, 4, 1) == pytest.approx(expected)
    assert bezier.cubicBezierPoint(t, -2, 2, 3, 1) == pytest.approx(expected)


@pytest.mark.parametrize('count', [1, 2, 3])
def test_path_segments_and_ownership(count):
    controls = list(range(3*count+1))
    path = bezier.BezierPath()
    path.setControlPoints(controls)
    assert path.getControlPoints() is controls
    assert path.curveCount == count and isinstance(path.curveCount, int)
    for i in range(count):
        for t in [0, 0.25, 0.5, 1]:
            assert path.calculateBezerPoint(i,t) == pytest.approx(3*i+3*t)
    assert controls == list(range(3*count+1))


def test_legacy_compatibility_and_integer_iteration(monkeypatch):
    from gem.experimental import bezier as legacy
    assert legacy.cubicBezierPoint is bezier.cubicBezierPoint
    assert legacy.quadraticBezierPoint is bezier.quadraticBezierPoint
    assert issubclass(legacy.BezierPath, bezier.BezierPath)
    path = legacy.BezierPath()
    path.setControlPoints(list(range(7)))
    # Isolate Python 3 iteration from unresolved adaptive sampling (E04).
    monkeypatch.setattr(path, 'findDrawingPoints', lambda i: [3*i, 3*i+3])
    assert path.getDrawingPoints() == [[0,3], [6]]
    assert path.calculateBezerPoint(1,0.5) == pytest.approx(4.5)
    other = legacy.BezierPath()
    other.interpolate([0,3], 1/3)
    assert other.controlPoints == [0,1,2,3]
    assert other.curveCount == 1 and isinstance(other.curveCount,int)
