import pytest
from gem import plane, ray, matrix, quaternion, vector


def V(*xs):
    return vector.Vector(len(xs),list(xs))


def test_plane_from_coefficients():
    p = plane.Plane()
    p.fromCoeffs(0,2,0,-4)
    assert [p.a,p.b,p.c,p.d] == [0,2,0,-4]
    assert p.normal.vector == [0,2,0]


def test_plane_normalization_offset():
    assert plane.normalize([0,2,0,-4]) == pytest.approx([0,1,0,-2])


def test_plane_normalization_normal_field():
    p = plane.Plane()
    p.a,p.b,p.c,p.d = 0,2,0,-4
    p.normal = V(0,2,0)
    assert p.normalize().normal.vector == [0,1,0]


def test_plane_from_points_incidence():
    p = plane.Plane()
    p.fromPoints(V(0,2,0),V(0,2,1),V(1,2,0))
    assert p.dot(V(0,2,0,1)) == pytest.approx(0)


def test_best_fit_normal_closed_polygon():
    vertices = [V(0,0,0),V(1,0,0),V(1,1,0),V(0,1,0)]
    assert abs(plane.Plane().bestFitNormal(vertices).vector[2]) == 1


def test_plane_manual_coefficients():
    p = plane.Plane()
    p.a,p.b,p.c,p.d = 0,1,0,-2
    p.normal = V(0,1,0)
    assert p.dot(V(0,2,0,1)) == 0
    assert p.point_location(p,[0,3,0]) == 1
    assert p.point_location(p,[0,1,0]) == -1
    assert p.point_location(p,[0,2,0]) == 0
    assert p.flip().flip().d == p.d
    assert p.clone().normal is not p.normal
    assert p.bestFitD([V(0,2,0),V(0,2,1)],p.normal) == 2


def test_duplicate_independence_and_distance():
    r = ray.Ray(V(1,2,3),V(0,0,5))
    copy = r.duplicate()
    assert copy.distance == r.distance
    copy.start.vector[0] = 99
    assert r.start.vector[0] == 1


def test_ray_quaternion_rotation():
    r = ray.Ray(V(1,0,0),V(1,0,0))
    r.rotateUsingQuaternion(quaternion.quat_from_axis_angle(V(0,0,1),90))
    assert r.start.vector == pytest.approx([0,1,0],abs=1e-14)
    assert r.dir.vector == pytest.approx([0,1,0],abs=1e-14)


def test_ray_translation_moves_origin():
    r = ray.Ray(V(1,2,3),V(0,0,1))
    r.translate(matrix.Matrix(4).translate(V(2,3,4)))
    assert r.start.vector == [3,5,7]
    assert r.dir.vector == [0,0,1]


def test_ray_matrix_rotation_3d():
    r = ray.Ray(V(1,0,0),V(2,0,0))
    assert r.distance == 2
    r.roateUsingMatrix(matrix.Matrix(3).rotate(V(0,0,1),90))
    assert r.start.vector == pytest.approx([0,1,0],abs=1e-14)
    assert r.dir.vector == pytest.approx([0,1,0],abs=1e-14)
