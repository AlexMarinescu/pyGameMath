"""Execute with an isolated interpreter containing an installed gem wheel."""
import math
from pathlib import Path
import sys

from gem import bezier, legendre
from gem.vector import Vector
from gem.experimental.bezier import cubicBezierPoint
from gem.experimental.legendre import Legendre


def main():
    source_root = Path(__file__).resolve().parents[1]
    for module in (bezier, legendre):
        assert source_root not in Path(module.__file__).resolve().parents
        assert Path(sys.prefix) in Path(module.__file__).resolve().parents
    assert cubicBezierPoint is bezier.cubicBezierPoint
    assert Legendre is legendre.Legendre
    comparisons = 0
    for degree, t, expected in [(2,1e-200,1e-100), (3,1e-150,1e-150)]:
        function = bezier.quadraticBezierPoint if degree == 2 else bezier.cubicBezierPoint
        values = [0.]*degree+[1e300]
        assert math.isclose(function(t,*values),expected,rel_tol=3e-14,abs_tol=0)
        comparisons += 1
        for dimension in (2,3):
            rows = [[p*(-1)**j for j in range(dimension)] for p in values]
            controls = [Vector(dimension,row) for row in rows]
            result = function(t,*controls)
            assert all(math.isclose(actual,expected*(-1)**j,rel_tol=3e-14,abs_tol=0)
                       for j,actual in enumerate(result.vector))
            assert all(result.vector is not row for row in rows)
            comparisons += 1
    for x in (-1e154,1e154,-10**154,10**154):
        p = Legendre(2,0,x); snapshot = p.__dict__.copy()
        assert math.isclose(p.run(),1.5e308,rel_tol=3e-15)
        assert p.__dict__ == snapshot
        comparisons += 1
    print('{} installed-wheel range regressions passed'.format(comparisons))


if __name__ == '__main__':
    main()
