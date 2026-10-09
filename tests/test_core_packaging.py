"""Core/compatibility boundaries, including a clean offline distribution build.

Building distributions requires setuptools and wheel in the build interpreter;
GEM_BUILD_PYTHON can select it independently of the pytest environment.
"""
import ast
import importlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('name,symbol', [
    ('bezier', 'BezierPath'), ('legendre', 'Legendre'),
    ('spherical_harmonics', 'rotate_coefficients'),
])
def test_canonical_core_import(name, symbol):
    module = importlib.import_module('gem.' + name)
    assert getattr(module, symbol).__module__ == 'gem.' + name


def test_retired_transport_import_is_absent():
    with pytest.raises(ModuleNotFoundError) as error:
        importlib.import_module('gem.experimental.sph_object')
    assert error.value.name == 'gem.experimental.sph_object'
    assert not (ROOT / 'gem/experimental/sph_object.py').exists()


def imports(path):
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            yield from (alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom):
            yield node.module or ''
            yield from ((node.module or '') + '.' + alias.name for alias in node.names)


def test_no_runtime_or_example_depends_on_experimental():
    for folder in [ROOT / 'gem', ROOT / 'examples']:
        for path in folder.rglob('*.py'):
            assert not any(name.startswith('gem.experimental') for name in imports(path)), path


def test_compatibility_modules_only_reexport_core():
    for path in (ROOT / 'gem/experimental').glob('*.py'):
        for node in ast.parse(path.read_text()).body:
            if isinstance(node, ast.Expr):
                assert isinstance(node.value, ast.Constant) and isinstance(node.value.value, str)
            elif isinstance(node, ast.ImportFrom):
                assert node.module in ('gem.bezier', 'gem.legendre', 'gem.spherical_harmonics')
            else:
                assert isinstance(node, ast.Assign)
                assert len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)
                assert node.targets[0].id == '__all__'
                assert isinstance(ast.literal_eval(node.value), list)


def test_clean_wheel_and_sdist_isolated_imports(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    shutil.copytree(ROOT / 'gem', source / 'gem',
                    ignore=shutil.ignore_patterns('__pycache__', '*.pyc', '*.egg-info'))
    for name in ('setup.py', 'setup.cfg', 'MANIFEST.in', 'README.rst', 'LICENSE'):
        shutil.copy2(ROOT / name, source / name)
    (source / 'docs').mkdir()
    for guide in ('EXPERIMENTAL_MIGRATION.md', 'VECTOR_VIEWPORT_CONTRACTS.md'):
        shutil.copy2(ROOT / 'docs' / guide, source / 'docs')
    build_python = os.environ.get('GEM_BUILD_PYTHON', sys._base_executable)
    env = os.environ.copy()
    env.pop('PYTHONPATH', None)
    subprocess.run([build_python, 'setup.py', 'sdist', 'bdist_wheel'],
                   cwd=source, env=env, check=True, capture_output=True, text=True)
    wheel, = (source / 'dist').glob('*.whl')
    sdist, = (source / 'dist').glob('*.tar.gz')
    required = ['gem/bezier.py', 'gem/legendre.py', 'gem/spherical_harmonics.py',
                'gem/experimental/bezier.py', 'gem/experimental/legendre.py',
                'gem/experimental/sph.py', 'gem/experimental/sph_sample.py',
                'gem/experimental/sph_irradiance_map.py',
                'gem/experimental/_bezier_legacy.py', 'gem/experimental/__init__.py']
    with zipfile.ZipFile(wheel) as archive:
        names = archive.namelist()
        assert all(name in names for name in required)
        assert not any('sph_object' in name for name in names)
    with tarfile.open(sdist) as archive:
        names = archive.getnames()
        assert all(any(name.endswith('/' + item) for name in names) for item in required)
        assert not any('sph_object' in name for name in names)
        assert any(name.endswith('/docs/EXPERIMENTAL_MIGRATION.md') for name in names)
        assert any(name.endswith('/docs/VECTOR_VIEWPORT_CONTRACTS.md') for name in names)
    # Install using the existing tooling; no index, resolution or dependency installs.
    target = tmp_path / 'installed'
    subprocess.run([sys.executable, '-m', 'pip', 'install', '--no-index', '--no-deps',
                    '--target', str(target), str(wheel)], cwd=tmp_path, env=env,
                   check=True, capture_output=True, text=True)
    # -I ignores source/PYTHONPATH. Supply only installed package plus existing
    # six's directory (the established dependency) to the isolated interpreter.
    import six
    code = '''
import sys, importlib, math
sys.path[:0] = [sys.argv[1], sys.argv[2]]
from gem import bezier, legendre, spherical_harmonics as sh
from gem.experimental import bezier as b, legendre as l, sph, sph_sample, sph_irradiance_map
from gem.experimental._bezier_legacy import BezierPath
assert b.BezierPath is BezierPath is bezier.BezierPath
assert l.Legendre is sph.Legendre is legendre.Legendre
assert sph.SPH is sh.SPH and sph.K is sh.K and sph.Factorial is sh.Factorial
assert sph_sample.SPHSample is sh.SPHSample and sph_sample.GenerateSamples is sh.GenerateSamples
assert sph_irradiance_map.SPH_IrradianceMapCoeff is sh.SPH_IrradianceMapCoeff
assert bezier.quadraticBezierPoint(.5, 0, 2, 4) == 2
assert abs(legendre.Legendre(3, 0, .2).run() + .28) < 1e-15
assert abs(sh.SPH(0, 0, 0, 0) - 1/math.sqrt(4*math.pi)) < 1e-15
from gem.vector import Vector, clamp
from gem.common import getViewPort
assert Vector(0) == Vector(0) and not (Vector(0) != Vector(0))
assert Vector(2) != Vector(3) and Vector(3) != Vector(2)
values = [-2, 2, 10]
clamped = clamp(3, values, [0]*3, [5]*3)
assert clamped.vector == [0, 2, 5] and values == [-2, 2, 10]
assert clamped.vector is not values
assert getViewPort(Vector(2, [3, 4]), 100, 200) == [83, 184, 100, 200]
for module in (bezier, legendre, sh, b, l, sph, sph_sample, sph_irradiance_map):
    assert module.__file__.startswith(sys.argv[1])
try:
    importlib.import_module('gem.experimental.sph_object')
except ModuleNotFoundError as error:
    assert error.name == 'gem.experimental.sph_object'
else:
    raise AssertionError('Retired transport shipped')
print('isolated wheel imports and known answers passed')
'''
    result = subprocess.run([sys.executable, '-I', '-c', code, str(target),
                             str(Path(six.__file__).parent)], cwd=tmp_path,
                            env=env, check=True, capture_output=True, text=True)
    assert 'known answers passed' in result.stdout
