"""Audit-only documentation/distribution checks; no installed runtime code.

Examples (from outside the checkout, replace the absolute checkout path):
  python -I /path/to/audit/phase4g5_verify.py smoke --output /tmp/wheel.json
  python /path/to/audit/phase4g5_verify.py docs --output /tmp/docs.json

The historical docs checker freezes a pre-repair checkpoint. This adapter
explicitly supplies this audit's base for its unchanged preservation checks.
"""
import argparse
from fractions import Fraction
import importlib
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fe453893b2c3a68b6619206cc3214036bf8c24c1'


def load_tool(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def documentation(package_root=None):
    checker = load_tool('audit_docs_checker', ROOT/'tools/check_architecture_docs.py')
    original_base = checker.BASE
    checker.BASE = BASE
    result = checker.check(examples=True, getting_started=True, package_root=package_root,
                           api_reference=True, tutorials=True, showcase=True, website=True)
    result['legacy_tool_checkpoint'] = original_base
    result['audit_checkpoint_explicitly_supplied'] = BASE
    result['inventory'] = checker.declarations()
    return result


def close(actual, expected, tolerance=3e-14):
    assert len(actual) == len(expected)
    for a,b in zip(actual,expected):
        assert math.isclose(a,b,rel_tol=tolerance,abs_tol=tolerance), (actual,expected)


def smoke():
    import gem
    package_root = Path(gem.__file__).resolve().parent.parent
    assert package_root.is_relative_to(Path(sys.prefix).resolve())
    assert not package_root.is_relative_to(ROOT)
    from gem import bezier, legendre, matrix, plane, quaternion, ray, spherical_harmonics as sh
    from gem.vector import Vector, refract
    modules = sorted((ROOT/'gem').glob('*.py'))
    imported = {}
    for path in modules:
        module = importlib.import_module('gem.'+path.stem) if path.stem != '__init__' else gem
        assert Path(module.__file__).resolve().is_relative_to(package_root)
        assert Path(module.__file__).read_bytes() == path.read_bytes()
        imported[module.__name__] = str(Path(module.__file__).resolve())
    aliases = {'bezier': ('bezier', ['BezierPath','cubicBezierPoint','quadraticBezierPoint']),
               'legendre': ('legendre', ['Legendre']),
               'sph': ('spherical_harmonics', ['Factorial','K','SPH','Legendre']),
               'sph_sample': ('spherical_harmonics', ['SPHSample','GenerateSamples']),
               'sph_irradiance_map': ('spherical_harmonics', ['SPH_IrradianceMapCoeff'])}
    for old,(new,names) in aliases.items():
        shim = importlib.import_module('gem.experimental.'+old)
        core = importlib.import_module('gem.'+new)
        for name in names:
            assert getattr(shim,name) is getattr(core,name)
    try:
        importlib.import_module('gem.experimental.sph_object')
    except ModuleNotFoundError as error:
        assert error.name == 'gem.experimental.sph_object'
    else:
        raise AssertionError('retired transport was installed')
    checks = []
    def passed(name):
        checks.append(name)
    q = quaternion.quat_from_axis_angle([0.,0.,1.],90.)
    close(quaternion.quat_rotate_vector(q,Vector(3,[1.,0.,0.])).vector,[0.,1.,0.])
    close((q.toMatrix()*Vector(4,[1.,0.,0.,0.])).vector,[0.,1.,0.,0.])
    close(quaternion.quat_from_matrix(q.toMatrix()).data,q.data)
    passed('axis-angle / Hamilton rotation / matrix / matrix-to-quaternion')
    point = Vector(3,[1.,2.,3.]); direction = Vector(3,[1.,0.,0.])
    r = ray.Ray(point,direction)
    r.distance = 7.
    r.rotateUsingQuaternion(q)
    close(r.start.vector,[-2.,1.,3.]); close(r.dir.vector,[0.,1.,0.])
    r.translate(matrix.Matrix(4).translate(Vector(3,[2.,-3.,4.])))
    close(r.start.vector,[0.,-2.,7.]); assert r.distance == 7.
    clone = r.duplicate()
    assert clone.start.vector is not r.start.vector and clone.end.vector is not r.end.vector
    assert point.vector == [1.,2.,3.]
    passed('ray rotation / translation / deep duplicate / constructor ownership')
    p = plane.Plane(); p.fromPoints(Vector(3,[0.,0.,2.]),Vector(3,[1.,0.,2.]),Vector(3,[0.,1.,2.]))
    assert (p.a,p.b,p.c,p.d) == (0.,0.,1.,-2.)
    assert p.dot(Vector(4,[2.,3.,2.,1.])) == 0.
    assert p.bestFitNormal([Vector(3,[1e16,1e16,2.]),Vector(3,[1e16+4,1e16,2.]),
                            Vector(3,[1e16+4,1e16+4,2.]),Vector(3,[1e16,1e16+4,2.])]).vector == [0.,0.,1.]
    incident = Vector(3,[1.,-1e-9,0.])
    assert refract(1.,incident,Vector(3,[0.,1.,0.])).vector == incident.vector
    close(refract(2/3.,Vector(3,[.6,-.8,0.]),Vector(3,[0.,1.,0.])).vector,
          [.4,-math.sqrt(.84),0.])
    passed('plane incidence / translated polygon normal / refraction')
    for size in (3,4):
        for scale in (1e-300,1e300):
            m = matrix.Matrix(size,[[scale*(i+1) if i==j else 0. for j in range(size)] for i in range(size)])
            inv = m.inverse()
            for i in range(size):
                for j in range(size):
                    expected = 1/scale/(i+1) if i==j else 0.
                    assert inv.matrix[i][j] == expected or math.isclose(inv.matrix[i][j],expected,rel_tol=3e-15)
                    import ctypes
                    assert inv.c_matrix[i][j] == ctypes.c_float(inv.matrix[i][j]).value
            for product in (m*inv,inv*m):
                for i,row in enumerate(product.matrix): close(row,[float(i==j) for j in range(size)])
    passed('scaled inverses / both multiplication orders / ctypes')
    assert bezier.quadraticBezierPoint(1e-200,0.,0.,1e300) == float(Fraction(1e-200)**2*Fraction(1e300))
    close(bezier.cubicBezierPoint(.5,Vector(2,[0.,0.]),Vector(2,[1.,2.]),
                                  Vector(2,[2.,2.]),Vector(2,[3.,0.])).vector,[1.5,1.5])
    path = bezier.BezierPath(); path.setControlPoints([Vector(2,[float(i),0.]) for i in range(7)])
    assert [[p.vector for p in segment] for segment in path.getDrawingPoints()] == [[[0.,0.],[3.,0.]],[[6.,0.]]]
    passed('scalar/vector Bezier / adaptive connected path')
    expected = float((3*Fraction(1e154)**2-1)/2)
    assert math.isclose(legendre.Legendre(2,0,1e154).run(),expected,rel_tol=3e-15)
    assert legendre.Legendre(2,1,.5).run() == -1.5*math.sqrt(.75)
    passed('Legendre finite extrapolation / Condon-Shortley phase')
    coefficients = [[0.]*3 for _ in range(9)]
    coefficients[0] = [math.sqrt(4*math.pi)]*3
    coefficients[3] = [-math.sqrt(4*math.pi/3),0.,0.]
    rotated = sh.rotate_coefficients(coefficients,q)
    close(sh.reconstruct(sh.convolve_diffuse(rotated),[0.,1.,0.]),[5*math.pi/3,math.pi,math.pi])
    assert math.isclose(sh.SPH(1,1,1e-9,0.),-math.sqrt(3/(4*math.pi))*math.sin(1e-9),rel_tol=3e-15)
    passed('SH canonical basis / near pole / analytical rotation / diffuse convolution')
    u = math.ulp(0.)
    endpoint = quaternion.quat_from_axis_angle([1.,0.,0.],math.degrees(2*u))
    actual = quaternion.quat_rotate_vector(quaternion.quat_slerp(quaternion.Quaternion(),endpoint,1.),
                                          Vector(3,[0.,1.,0.])).vector
    expected = [0.,1.,float(2*Fraction(endpoint.data[0])*Fraction(endpoint.data[1]))]
    assert actual == [0.,1.,0.] and expected[2] == 2*u
    return {'package_root': str(package_root), 'core_imports': imported,
            'compatibility_modules': aliases, 'retired_transport_absent': True,
            'passed_checks': checks, 'mathematical_smoke_checks_passed': len(checks),
            'confirmed_finding': {'id': '4G5-A01', 'reproduced': True, 'actual': actual, 'expected': expected},
            'documentation': documentation(package_root),
            'dependencies': {name: importlib.metadata.version(name) for name in ('gem','six')}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=('smoke','docs'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = smoke() if args.mode == 'smoke' else documentation()
    result['environment'] = {'python': sys.version, 'platform': platform.platform(),
                             'machine': platform.machine(), 'executable': sys.executable}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    print('Audit verification recorded in '+str(args.output))


if __name__ == '__main__':
    main()
