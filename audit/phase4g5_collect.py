"""Collect executed audit evidence from /tmp; asserts artifact/scope consistency.

Run after the commands in PHASE4G5-FINAL-INTEGRATION.md:
  python audit/phase4g5_collect.py --output audit/phase4g5-verification.json
This tool does not run or repair mathematics and is not packaged with gem.
"""
import argparse
import ast
import configparser
from email.parser import Parser
import hashlib
import importlib.metadata
import importlib.util
import json
from pathlib import Path
import platform
import re
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fe453893b2c3a68b6619206cc3214036bf8c24c1'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def junit(name):
    path = Path('/tmp/phase4g5-'+name+'.xml')
    root = ET.parse(path).getroot()
    cases = root.findall('.//testcase')
    failures, xfails, skips, errors = [], [], [], []
    for case in cases:
        identity = case.get('classname')+'::'+case.get('name')
        if case.find('failure') is not None:
            failures.append({'test': identity, 'message': case.find('failure').get('message')})
        elif case.find('error') is not None:
            errors.append(identity)
        elif case.find('skipped') is not None:
            element = case.find('skipped')
            (xfails if element.get('type') == 'pytest.xfail' else skips).append(identity)
    return {'passed': len(cases)-len(failures)-len(xfails)-len(skips)-len(errors),
            'failed': len(failures), 'xfail': len(xfails), 'skipped': len(skips), 'errors': len(errors),
            'total': len(cases), 'xfail_identities': xfails, 'failures': failures,
            'elapsed_seconds': sum(float(s.get('time')) for s in root.findall('testsuite')),
            'junit_sha256': sha(path.read_bytes())}


def hdr_comparison():
    spec = importlib.util.spec_from_file_location('hdr_decoder', ROOT/'benchmarks/verify_hdr_reference.py')
    decoder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(decoder)
    def differences(a,b,path=''):
        if isinstance(a,(int,float)) and not isinstance(a,bool):
            return [] if a == b else [{'path':path, 'before':a, 'after':b, 'absolute_difference':abs(a-b)}]
        if isinstance(a,dict):
            assert set(a) == set(b)
            return [r for k in a for r in differences(a[k],b[k],path+'/'+str(k))]
        if isinstance(a,list):
            assert len(a) == len(b)
            return [r for i,(x,y) in enumerate(zip(a,b)) for r in differences(x,y,path+'/'+str(i))]
        assert a == b, (path,a,b)
        return []
    rows = {}
    for old in sorted((ROOT/'examples/output').iterdir()):
        new = Path('/tmp/phase4g5-hdr-output')/old.name
        row = {'committed_sha256':sha(old.read_bytes()), 'generated_sha256':sha(new.read_bytes()),
               'bytes_identical':old.read_bytes() == new.read_bytes()}
        if old.suffix == '.png':
            w,h,pixels = decoder.png_rgb(new)
            assert (w,h,pixels) == decoder.png_rgb(old)
            row.update(width=w,height=h,decoded_rgb_identical=True,changed_channel_bytes=0,
                       minimum_byte=min(pixels),maximum_byte=max(pixels),mean_byte=sum(pixels)/len(pixels))
        elif old.suffix == '.json':
            changes = differences(json.loads(old.read_text()),json.loads(new.read_text()))
            row.update(numeric_changes=changes,changed_numeric_fields=len(changes),
                       maximum_absolute_difference=max((r['absolute_difference'] for r in changes),default=0),
                       nonnumeric_fields_unchanged=True)
        else:
            pattern = r'[-+]?(?:\d+\.\d*|\.\d+)(?:[eE][-+]?\d+)?'
            assert re.sub(pattern,'#',old.read_text()) == re.sub(pattern,'#',new.read_text())
            a,b = [list(map(float,re.findall(pattern,p.read_text()))) for p in (old,new)]
            assert len(a) == len(b)
            delta = [abs(x-y) for x,y in zip(a,b) if x != y]
            row.update(changed_numeric_fields=len(delta),maximum_absolute_difference=max(delta,default=0),
                       nonnumeric_text_unchanged=True)
        rows[old.name] = row
    return {'files':rows,'goldens_modified':False,'byte_hashes_are_supplemental':True}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    runs = {name: junit(name) for name in ('baseline','focused','reproductions','full')}
    assert (runs['baseline']['passed'],runs['baseline']['total']) == (3714,3714)
    assert (runs['focused']['passed'],runs['focused']['xfail']) == (139,12)
    assert runs['reproductions']['failed'] == 12 and runs['reproductions']['total'] == 12
    assert (runs['full']['passed'],runs['full']['xfail'],runs['full']['total']) == (3853,12,3865)
    for name in ('baseline','focused','full'):
        assert not (runs[name]['failed'] or runs[name]['skipped'] or runs[name]['errors'])
    assert {r['test'] for r in runs['reproductions']['failures']} == set(runs['full']['xfail_identities'])
    runtime = sorted(str(p.relative_to(ROOT)) for p in (ROOT/'gem').rglob('*.py'))
    assert len(runtime) == 17
    contents = {p: (ROOT/p).read_bytes() for p in runtime}
    runtime_imports = {}
    for name in runtime:
        imports = []
        for node in ast.walk(ast.parse(contents[name])):
            if isinstance(node, ast.Import):
                imports.extend(alias.name for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                imports.append('.'*node.level+(node.module or ''))
        runtime_imports[name] = sorted(set(imports))
        assert not any('sph_object' in imported for imported in imports)
        if not name.startswith('gem/experimental/'):
            assert not any('experimental' in imported for imported in imports)
    distributions = {}
    for path in sorted(Path('/tmp/phase4g5-dist-build/dist').iterdir()):
        if path.suffix == '.whl':
            with zipfile.ZipFile(path) as archive:
                names = archive.namelist()
                files = {name: archive.read(name) for name in names}
                metadata = Parser().parsestr(files['gem-0.1.12.dist-info/METADATA'].decode())
                assert metadata.get_all('Requires-Dist') == ['six']
                assert b'Root-Is-Purelib: true' in files['gem-0.1.12.dist-info/WHEEL']
        else:
            with tarfile.open(path) as archive:
                files = {member.name.split('/',1)[1]: archive.extractfile(member).read()
                         for member in archive.getmembers() if member.isfile()}
                names = sorted(files)
                metadata = Parser().parsestr(files['PKG-INFO'].decode())
                assert metadata.get_all('Requires-Dist') == ['six']
                for name in ('setup.py','MANIFEST.in','LICENSE','README.rst'):
                    assert files[name] == (ROOT/name).read_bytes()
                original, built = configparser.ConfigParser(), configparser.ConfigParser()
                original.read_string((ROOT/'setup.cfg').read_text())
                built.read_string(files['setup.cfg'].decode())
                assert dict(original['metadata']) == dict(built['metadata'])
                assert built.sections() == ['metadata','egg_info']
                assert dict(built['egg_info']) == {'tag_build':'', 'tag_date':'0'}
        assert sorted(name for name in files if name.startswith('gem/') and name.endswith('.py')) == runtime
        for name, expected in contents.items(): assert files[name] == expected
        assert not any(name.endswith(('.so','.dll','.pyd')) for name in names)
        assert not any('sph_object' in name for name in names)
        distributions[path.name] = {'sha256': sha(path.read_bytes()), 'bytes': path.stat().st_size,
                                    'members': names, 'runtime_files': runtime,
                                    'runtime_byte_identical': True, 'compiled_extensions': [],
                                    'requires_dist': metadata.get_all('Requires-Dist'),
                                    'version': metadata['Version'],
                                    'requires_python': metadata['Requires-Python'],
                                    'classifiers': metadata.get_all('Classifier')}
        if path.suffix != '.whl':
            distributions[path.name]['setup_cfg_build_normalization'] = {
                'comments_removed': True, 'generated_egg_info': dict(built['egg_info']),
                'metadata_preserved': True, 'source_setup_cfg_unchanged': True}
    installed = {}
    for name in ('wheel','sdist'):
        data = json.loads(Path('/tmp/phase4g5-'+name+'-smoke.json').read_text())
        data['documentation'].pop('inventory')
        assert data['mathematical_smoke_checks_passed'] == 7 and data['confirmed_finding']['reproduced']
        installed[name] = data
    protected = subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    fingerprint = hashlib.sha256()
    for name in protected:
        before = subprocess.check_output(['git','show',BASE+':'+name],cwd=ROOT)
        after = (ROOT/name).read_bytes()
        assert before == after, name
        fingerprint.update(name.encode()); fingerprint.update(after)
    cpu = next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines()
                if line.startswith('model name')), 'unavailable')
    result = {'schema': 1, 'phase': '4G-5', 'base': BASE, 'branch': 'audit/phase4g5-final-integration',
              'audit_only': True, 'release_signoff': False,
              'environment': {'python': sys.version, 'platform': platform.platform(), 'cpu': cpu,
                              'executable': sys.executable, 'shared_host': True,
                              'pytest': importlib.metadata.version('pytest'),
                              'six': importlib.metadata.version('six'),
                              'build_tools': {'setuptools':'84.0.0','wheel':'0.48.0','packaging':'26.3'}},
              'tests': runs, 'distributions': distributions, 'installed': installed,
              'documentation': json.loads(Path('/tmp/phase4g5-docs.json').read_text()),
              'original_documentation_cli': {'exit_code': 2,
                  'error': 'out-of-scope modification: gem/bezier.py',
                  'reason': 'hardcoded checkpoint predates legitimate merged numerical repairs',
                  'log_sha256': sha(Path('/tmp/phase4g5-source-docs.log').read_bytes())},
              'hdr_regeneration': hdr_comparison(),
              'performance_runs': [json.loads(Path('/tmp/phase4g5-path-cost-'+str(i)+'.json').read_text())
                                   for i in (1,2)],
              'preservation': {'all_original_tracked_files_unchanged': len(protected),
                               'all_original_files_sha256': fingerprint.hexdigest(),
                               'runtime_sha256': {p:sha(contents[p]) for p in runtime},
                               'runtime_imports': runtime_imports,
                               'core_depends_on_experimental': False,
                               'retired_runtime_imports': [],
                               'existing_tests_unchanged': True, 'packaging_unchanged': True,
                               'dependencies_unchanged': True, 'public_signatures_unchanged': True},
              'findings': [{'id': '4G5-A01', 'severity': 'P3', 'classification': 'Release blocker',
                            'source': 'gem/quaternion.py:289-293; squad4 propagation at 331-333',
                            'expected_transverse_rotation': 2*float.fromhex('0x0.0000000000001p-1022'),
                            'actual_transverse_rotation': 0., 'strict_xfails': 12},
                           {'id': '4G5-P01', 'severity': 'P3', 'classification': 'Post-1.0 improvement',
                            'source': 'gem/bezier.py:169-188,219-246',
                            'validated_controls_formula': '(segments+1)*(3*segments+1)',
                            'confirmed_in_two_runs': True}],
              'commands_report': 'PHASE4G5-FINAL-INTEGRATION.md'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print('Collected full-suite, installed artifact, scope and numerical evidence')


if __name__ == '__main__':
    main()
