"""Inspect release archives and optional isolated installed imports; never upload."""
import argparse
import hashlib
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import tarfile
import zipfile
from email.parser import BytesParser

ROOT = Path(__file__).resolve().parents[1]


def digest(data):
    return hashlib.sha256(data).hexdigest()


def inspect_archives(wheel, sdist):
    runtime = {str(p.relative_to(ROOT)): p.read_bytes() for p in (ROOT/'gem').rglob('*.py')}
    with zipfile.ZipFile(wheel) as archive:
        files = {name: archive.read(name) for name in archive.namelist() if not name.endswith('/')}
    with tarfile.open(sdist) as archive:
        source = {member.name.split('/',1)[1]: archive.extractfile(member).read()
                  for member in archive.getmembers() if member.isfile()}
    for name, content in runtime.items():
        assert files[name] == source[name] == content, name
    assert {n for n in files if n.startswith('gem/')} == set(runtime)
    for names in (files, source):
        assert not any('__pycache__' in n or n.endswith(('.pyc','.pyo')) for n in names)
        assert not any(n.startswith(('build/', 'dist/', 'site/', '.venv', '.git/')) for n in names)
        assert not any('sph_object' in n for n in names)
    metadata_name, = [n for n in files if n.endswith('.dist-info/METADATA')]
    metadata = BytesParser().parsebytes(files[metadata_name])
    source_metadata = BytesParser().parsebytes(source['PKG-INFO'])
    version = {}
    exec((ROOT/'gem/_version.py').read_text(), version)
    for message in (metadata, source_metadata):
        assert message['Name'] == 'gem'
        assert message['Version'] == version['__version__'] == '1.0.0'
        assert message['Requires-Python'] == '>=3.10'
        assert message.get_all('Requires-Dist') == ['six']
        assert message['License-Expression'] == 'BSD-2-Clause'
        assert message.get_all('License-File') == ['LICENSE']
        assert 'explosiveduck' not in str(message)
        assert set(message.get_all('Project-URL')) == {
            'Homepage, https://github.com/AlexMarinescu/pyGameMath',
            'Documentation, https://alexmarinescu.github.io/pyGameMath/',
            'Source, https://github.com/AlexMarinescu/pyGameMath',
            'Issues, https://github.com/AlexMarinescu/pyGameMath/issues'}
        assert 'Programming Language :: Python :: 2.7' not in message.get_all('Classifier')
    assert not any(n.endswith('entry_points.txt') for n in files)
    assert any(n.endswith('/licenses/LICENSE') and files[n] == (ROOT/'LICENSE').read_bytes() for n in files)
    for name in ('LICENSE', 'README.md', 'README.rst', 'pyproject.toml',
                 'docs/development/packaging.md', 'docs/architecture/conventions.md',
                 'examples/hdr_sh/reference.py'):
        assert source[name] == (ROOT/name).read_bytes(), name
    return {'metadata': dict(metadata.items()), 'runtime_files': len(runtime),
            'runtime_sha256': {name:digest(data) for name,data in runtime.items()},
            'wheel': {'name':wheel.name, 'sha256':digest(wheel.read_bytes()), 'files':sorted(files)},
            'sdist': {'name':sdist.name, 'sha256':digest(sdist.read_bytes()), 'files':sorted(source)}}


def inspect_installed(package_root):
    sys.path.insert(0, str(package_root.resolve()))
    import gem
    assert Path(gem.__file__).resolve().is_relative_to(package_root.resolve())
    assert gem.__version__ == importlib.metadata.version('gem') == '1.0.0'
    runtime = list((ROOT/'gem').rglob('*.py'))
    for path in runtime:
        assert (package_root/path.relative_to(ROOT)).read_bytes() == path.read_bytes(), path
    for name in ('bezier', 'legendre', 'spherical_harmonics', 'common','vector','matrix','quaternion','ray','plane'):
        imported = importlib.import_module('gem.'+name)
        assert Path(imported.__file__).resolve().is_relative_to(package_root.resolve())
    from gem import bezier,legendre,spherical_harmonics as sh
    from gem.experimental import bezier as legacy_b,legendre as legacy_l,sph
    assert legacy_b.BezierPath is bezier.BezierPath
    assert legacy_l.Legendre is legendre.Legendre and sph.SPH is sh.SPH
    from gem.vector import Vector
    from gem.matrix import Matrix
    from gem.quaternion import Quaternion, squad4
    import math
    assert Vector(3,[1,2,3])+Vector(3,[4,5,6]) == Vector(3,[5,7,9])
    assert Vector(3,[1e300,-2e300,2e300]).normalize().vector == [1/3,-2/3,2/3]
    for size in (3,4):
        obj = Matrix(size, [[(2. if i==j else 1. if j==i+1 else 0.)*1e-300 for j in range(size)] for i in range(size)])
        result = obj.inverse()
        for i in range(size):
            for j in range(size):
                expected = ((-1)**(j-i)/2**(j-i+1))/1e-300 if j>=i else 0.
                assert math.isclose(result.matrix[i][j],expected,rel_tol=3e-14,abs_tol=0)
    tiny = float.fromhex('0x0.0000000000001p-1022')
    a,b = Quaternion([1.,0.,0.,0.]),Quaternion([1.,tiny,0.,0.])
    assert a.slerp(b,1).data == b.data
    assert squad4(a,b,a,b,1).data == b.data
    assert bezier.quadraticBezierPoint(.5,0.,2.,4.) == 2.
    assert abs(legendre.Legendre(3,0,.2).run()+.28)<1e-15
    assert sh.SPH(1,1,1e-9,0) != 0
    return {'python':platform.python_version(), 'implementation':platform.python_implementation(),
            'platform':platform.platform(), 'package_root':str(package_root.resolve()),
            'version':gem.__version__, 'runtime_files_matched':len(runtime), 'smoke':'passed'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheel',type=Path); parser.add_argument('--sdist',type=Path)
    parser.add_argument('--installed-root',type=Path); parser.add_argument('--output',type=Path)
    args=parser.parse_args(); report={}
    if args.wheel and args.sdist: report['artifacts']=inspect_archives(args.wheel,args.sdist)
    if args.installed_root: report['installed']=inspect_installed(args.installed_root)
    if not report: parser.error('supply wheel/sdist or installed-root')
    text=json.dumps(report,indent=2)+'\n'
    if args.output: args.output.write_text(text)
    print(text)


if __name__=='__main__': main()
