"""Collect executed repair evidence; run after the commands in the repair report."""
import argparse
import ast
import configparser
from email.parser import Parser
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import tarfile
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = '4e48dc68cc23d47abd7501897c6dfbb0c05156c3'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def before(name):
    return subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT)


def junit(name):
    path = Path('/tmp/phase4g5r-'+name+'.xml')
    root = ET.parse(path).getroot()
    groups = {key: [] for key in ('passed', 'failed', 'xfail', 'skipped', 'errors')}
    for case in root.findall('.//testcase'):
        identity = case.get('classname')+'::'+case.get('name')
        key = 'passed'
        if case.find('failure') is not None:
            key = 'failed'
        elif case.find('error') is not None:
            key = 'errors'
        elif case.find('skipped') is not None:
            key = 'xfail' if case.find('skipped').get('type') == 'pytest.xfail' else 'skipped'
        groups[key].append(identity)
    return {**{key: len(value) for key, value in groups.items()},
            'total': sum(map(len, groups.values())), 'identities': groups,
            'elapsed_seconds': sum(float(s.get('time')) for s in root.findall('testsuite')),
            'junit_sha256': sha(path.read_bytes())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    runs = {name: junit(name) for name in
            ('baseline', 'reproductions', 'original-passing', 'focused', 'full')}
    assert (runs['baseline']['passed'], runs['baseline']['xfail']) == (3853, 12)
    assert runs['reproductions']['failed'] == runs['reproductions']['total'] == 12
    assert runs['original-passing']['passed'] == runs['original-passing']['total'] == 12
    assert runs['focused']['passed'] == runs['focused']['total'] == 1120
    assert runs['full']['passed'] == runs['full']['total'] == 3947
    for name in ('original-passing', 'focused', 'full'):
        assert not any(runs[name][key] for key in ('failed', 'xfail', 'skipped', 'errors'))
    failed_ids = set(runs['reproductions']['identities']['failed'])
    assert failed_ids == set(runs['baseline']['identities']['xfail'])
    assert failed_ids == set(runs['original-passing']['identities']['passed'])
    assert failed_ids <= set(runs['full']['identities']['passed'])
    marker = b"@pytest.mark.defect('4G5-A01: subnormal SLERP separation loses the endpoint rotation')\n"
    original_tests = 'tests/test_final_integration_audit.py'
    assert before(original_tests).count(marker) == 1
    assert before(original_tests).replace(marker, b'') == (ROOT/original_tests).read_bytes()
    changed = {'gem/quaternion.py', original_tests, 'docs/api/quaternion.md'}
    protected = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE],
                                        cwd=ROOT, text=True).splitlines()
    fingerprint = hashlib.sha256()
    for name in protected:
        if name in changed:
            continue
        content = (ROOT/name).read_bytes()
        assert before(name) == content, name
        fingerprint.update(name.encode()); fingerprint.update(content)
    old_ast = ast.parse(before('gem/quaternion.py'))
    new_ast = ast.parse((ROOT/'gem/quaternion.py').read_bytes())
    def without_slerp(tree):
        return [ast.dump(node) for node in tree.body if getattr(node, 'name', '') != 'quat_slerp']
    assert without_slerp(old_ast) == without_slerp(new_ast)
    old_function = next(n for n in old_ast.body if getattr(n, 'name', '') == 'quat_slerp')
    new_function = next(n for n in new_ast.body if getattr(n, 'name', '') == 'quat_slerp')
    assert ast.dump(old_function.args) == ast.dump(new_function.args)
    runtime = sorted(str(p.relative_to(ROOT)) for p in (ROOT/'gem').rglob('*.py'))
    assert len(runtime) == 17
    contents = {p: (ROOT/p).read_bytes() for p in runtime}
    distributions = {}
    for path in sorted(Path('/tmp/phase4g5r-final-dist-build/dist').iterdir()):
        if path.suffix == '.whl':
            with zipfile.ZipFile(path) as archive:
                files = {name: archive.read(name) for name in archive.namelist()}
            metadata = Parser().parsestr(files['gem-0.1.12.dist-info/METADATA'].decode())
            assert b'Root-Is-Purelib: true' in files['gem-0.1.12.dist-info/WHEEL']
        else:
            with tarfile.open(path) as archive:
                files = {m.name.split('/', 1)[1]: archive.extractfile(m).read()
                         for m in archive.getmembers() if m.isfile()}
            metadata = Parser().parsestr(files['PKG-INFO'].decode())
            for name in ('setup.py', 'MANIFEST.in', 'LICENSE', 'README.rst'):
                assert files[name] == (ROOT/name).read_bytes()
            source, built = configparser.ConfigParser(), configparser.ConfigParser()
            source.read_string((ROOT/'setup.cfg').read_text())
            built.read_string(files['setup.cfg'].decode())
            assert dict(source['metadata']) == dict(built['metadata'])
            assert built.sections() == ['metadata', 'egg_info']
            assert dict(built['egg_info']) == {'tag_build': '', 'tag_date': '0'}
        assert metadata.get_all('Requires-Dist') == ['six']
        assert sorted(n for n in files if n.startswith('gem/') and n.endswith('.py')) == runtime
        assert all(files[name] == value for name, value in contents.items())
        assert not any(n.endswith(('.so', '.dll', '.pyd')) or 'sph_object' in n for n in files)
        distributions[path.name] = {'sha256': sha(path.read_bytes()), 'bytes': path.stat().st_size,
            'runtime_files': runtime, 'runtime_byte_identical': True, 'compiled_extensions': [],
            'requires_dist': metadata.get_all('Requires-Dist'), 'version': metadata['Version'],
            'setup_cfg_build_normalization': path.suffix != '.whl'}
    assert len(distributions) == 2
    smoke = {}
    for kind in ('source', 'wheel', 'sdist'):
        data = json.loads(Path('/tmp/phase4g5r-'+kind+'-smoke.json').read_text())
        assert data['original_regressions_passed'] == 12
        assert data['additional_subnormal_checks_passed'] == 6
        assert data['rotation_matrix_checks_passed'] == 2
        assert data['examples_executed'] == 42
        assert data['source_mode'] == (kind == 'source')
        if kind != 'source':
            assert 'No broken requirements found.' in Path('/tmp/phase4g5r-'+kind+'-pip-check.log').read_text()
        smoke[kind] = data
    performance = []
    for suffix in ('', '-repeat', '-third'):
        name = 'audit/phase4g5r-performance'+suffix+'.json'
        data = json.loads((ROOT/name).read_text())
        assert data['source_sha256']['after'] == sha(contents['gem/quaternion.py'])
        assert data['source_sha256']['before'] == sha(before('gem/quaternion.py'))
        assert data['base'] == BASE
        assert data['compatibility']['changed_endpoint_calls'] == 0
        performance.append({'file': name, 'sha256': sha((ROOT/name).read_bytes()),
                            'compatibility': data['compatibility'], 'method': data['method']})
    new_cases = [case for case in runs['full']['identities']['passed']
                 if case.startswith('tests.test_slerp_numerical_repair::')]
    assert len(new_cases) == 82
    for run in runs.values():
        passed = run['identities'].pop('passed')
        run['passed_identities_sha256'] = sha('\n'.join(sorted(passed)).encode())
        run['original_passing_regression_identities'] = sorted(failed_ids.intersection(passed))
    cpu = next((s.split(':', 1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines()
                if s.startswith('model name')), 'unavailable')
    result = {'schema': 1, 'phase': '4G-5R', 'base': BASE,
        'branch': 'repair/phase4g5-slerp-endpoints', 'finding': '4G5-A01',
        'environment': {'python': sys.version, 'platform': platform.platform(), 'cpu': cpu,
            'shared_host': True, 'executable': sys.executable,
            'pytest': importlib.metadata.version('pytest'), 'six': importlib.metadata.version('six'),
            'build_tools': {'setuptools': '84.0.0', 'wheel': '0.48.0', 'packaging': '26.3'},
            'other_python_versions_verified': False},
        'tests': runs, 'new_independent_cases': new_cases, 'installed_and_source': smoke,
        'distributions': distributions, 'performance_runs': performance,
        'preservation': {'original_regression_bodies_inputs_assertions_tolerances_unchanged': True,
            'removed_strict_xfail_markers': 1, 'removed_strict_xfail_cases': 12,
            'original_tracked_files_unchanged': len(protected)-len(changed),
            'unchanged_files_fingerprint_sha256': fingerprint.hexdigest(),
            'changed_original_files': sorted(changed), 'only_modified_runtime_function': 'quat_slerp',
            'other_runtime_files_unchanged': 16, 'public_signatures_unchanged': True,
            'packaging_dependencies_versions_unchanged': True, 'documentation_checker_unchanged': True,
            'bezier_implementation_unchanged': True,
            'runtime_sha256': {name: sha(content) for name, content in contents.items()}},
        'limitations': ['Binary64 product rounding remains; unrepresentable orientations are not recovered.',
            'Unit-input subnormal-angle limit is supported for t in [0,1]; no nonunit-domain expansion.',
            'Shared-host measurements are not universal performance guarantees.',
            'Documentation checks reuse unchanged signature/example helpers, not its historical scope gate.'],
        'commands_report': 'PHASE4G5R-SLERP-REPAIR.md'}
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print('Collected final tests, original assertions, distributions, examples and matched measurements')


if __name__ == '__main__':
    main()
