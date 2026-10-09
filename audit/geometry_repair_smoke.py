"""Run with an isolated installed interpreter, outside the source tree."""
from decimal import Decimal, localcontext
from fractions import Fraction
import importlib.metadata
import json
import math
from pathlib import Path
import sys

from gem import bezier, legendre, matrix, plane, quaternion, ray, spherical_harmonics, vector
from gem.experimental import bezier as old_bezier, legendre as old_legendre
from gem.experimental import sph as old_sh


def V(data):
    return vector.Vector(len(data), list(data))


def main():
    source = Path(__file__).resolve().parents[1]
    assert source not in [Path(p).resolve() for p in sys.path if p]
    for module in (bezier, legendre, matrix, plane, quaternion, ray, spherical_harmonics, vector):
        assert source not in Path(module.__file__).resolve().parents
    assert old_bezier.cubicBezierPoint is bezier.cubicBezierPoint
    assert old_legendre.Legendre is legendre.Legendre
    assert old_sh.SPH is spherical_harmonics.SPH
    cases = 0
    for size in (2, 3, 4):
        for component in (1e-9, 1e-100):
            values = [math.sqrt(1.-component*component), -component] + [0.]*(size-2)
            incident, normal = V(values), V([0., 1.] + [0.]*(size-2))
            storage = incident.vector
            result = vector.refract(1., incident, normal)
            assert result.vector == values and incident.vector is storage
            assert result.vector is not storage and result is not incident
            cases += 1
    for oblique in (False, True):
        for reverse in (False, True):
            for closed in (False, True):
                local = [[0., 0., 0.], [1., 0., float(oblique)], [0., 1., float(oblique)]]
                values = [[2.**52+x for x in point] for point in local]
                if reverse:
                    values.reverse()
                if closed:
                    values.append(values[0][:])
                with localcontext() as ctx:
                    ctx.prec = 110
                    length = Decimal(1+2*int(oblique)).sqrt()
                    sign = -1 if reverse else 1
                    expected = [float(Decimal(-int(oblique)*sign)/length)]*2 + [float(Decimal(sign)/length)]
                points = [V(point) for point in values]
                storages = [point.vector for point in points]
                result = plane.Plane().bestFitNormal(points)
                assert all(abs(a-b) <= 3e-15 for a,b in zip(result.vector, expected))
                assert all(p.vector is s and p.vector == v for p,s,v in zip(points, storages, values))
                cases += 1
    assert vector.refract(1., V([1., -5e-324]), V([0., 1.])).vector == [1., -5e-324]
    assert vector.refract(1.5, V([.9, -math.sqrt(.19), 0.]), V([0., 1., 0.])).vector == [0., 0., 0.]
    for points in ([], [V([0., 0., 0.])], [V([0., 0., 0.]), V([1., 1., 1.])]):
        try:
            plane.Plane().bestFitNormal(points)
        except ZeroDivisionError:
            pass
        else:
            raise AssertionError('degenerate polygon no longer raises ZeroDivisionError')
    # Existing supported core functionality and transitional imports remain usable.
    assert spherical_harmonics.SPH(1, 1, 1e-9, 0.) != 0.
    assert matrix.Matrix(3).inverse().matrix == matrix.Matrix(3).matrix
    assert quaternion.Quaternion([1., 0., 0., 0.]).pow(2).data == [1., 0., 0., 0.]
    assert Fraction(vector.refract(1., V([1., 0.]), V([0., 1.])).vector[0]) == 1
    print(json.dumps(dict(status='passed', original_regressions=cases,
                         gem_version=importlib.metadata.version('gem'),
                         vector_path=vector.__file__, plane_path=plane.__file__,
                         isolated=True, python=sys.version), indent=2))


if __name__ == '__main__':
    main()
