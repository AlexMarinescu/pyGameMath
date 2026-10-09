"""Collect executed-test evidence and check source/artifact preservation."""
import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile
import xml.etree.ElementTree as ET
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE = '495dc3d5b4d5f8bd4969b690b3a207a9f3339023'


def baseline(path):
    return subprocess.check_output(['git', 'show', BASE + ':' + str(path)], cwd=ROOT)


def junit(path):
    counts = Counter(passed=0, failed=0, xfailed=0, skipped=0, errors=0)
    cases = []
    for case in ET.parse(path).iter('testcase'):
        status = 'passed'
        if case.find('failure') is not None:
            status = 'failed'
        elif case.find('error') is not None:
            status = 'errors'
        elif case.find('skipped') is not None:
            skipped = case.find('skipped')
            status = 'xfailed' if skipped.get('type') == 'pytest.xfail' else 'skipped'
        counts[status] += 1
        cases.append(dict(name=case.get('name'), classname=case.get('classname'), status=status))
    return dict(counts), cases


def signatures(source):
    return [ast.dump(node.args) for node in ast.walk(ast.parse(source))
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evidence-dir', type=Path, default=Path('/tmp'))
    parser.add_argument('--dist-dir', type=Path, default=Path('/tmp/phase4g4r-dist-build/dist'))
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    evidence = args.evidence_dir
    baseline_counts, _ = junit(evidence/'phase4g4r-baseline.xml')
    original_counts, original_cases = junit(evidence/'phase4g4r-original.xml')
    full_counts, full_cases = junit(evidence/'phase4g4r-full.xml')
    assert original_counts == dict(passed=0, failed=14, xfailed=0, skipped=0, errors=0)
    assert baseline_counts == dict(passed=3540, failed=0, xfailed=14, skipped=0, errors=0)
    assert full_counts == dict(passed=3714, failed=0, xfailed=0, skipped=0, errors=0)
    ids = {(c['classname'], c['name']) for c in original_cases}
    repaired = [c for c in full_cases if (c['classname'], c['name']) in ids]
    assert len(repaired) == 14 and all(c['status'] == 'passed' for c in repaired)
    original_ast = ast.parse(baseline('tests/test_geometry_audit.py'))
    targets = {'test_equal_media_grazing_refraction_identity',
               'test_translated_polygon_retains_exact_area_normal'}
    removed = 0
    for node in ast.walk(original_ast):
        if isinstance(node, ast.FunctionDef) and node.name in targets:
            decorators = node.decorator_list
            node.decorator_list = [d for d in decorators if not (
                isinstance(d, ast.Call) and isinstance(d.func, ast.Attribute) and d.func.attr == 'defect')]
            removed += len(decorators)-len(node.decorator_list)
    assert removed == 2
    assert ast.dump(original_ast) == ast.dump(ast.parse((ROOT/'tests/test_geometry_audit.py').read_bytes()))
    runtime = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE, 'gem'], cwd=ROOT).decode().splitlines()
    runtime = [p for p in runtime if p.endswith('.py')]
    repaired_paths = {'gem/vector.py', 'gem/plane.py'}
    untouched = [p for p in runtime if p not in repaired_paths]
    assert all(baseline(p) == (ROOT/p).read_bytes() for p in untouched)
    assert all(signatures(baseline(p)) == signatures((ROOT/p).read_bytes()) for p in repaired_paths)
    # Check that only the two intended function bodies changed, not adjacent math.
    for path, name in [('gem/vector.py', 'refract'), ('gem/plane.py', 'bestFitNormal')]:
        old, new = ast.parse(baseline(path)), ast.parse((ROOT/path).read_bytes())
        for tree in (old, new):
            matches = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name]
            assert len(matches) == 1
            matches[0].body = [ast.Pass()]
        assert ast.dump(old) == ast.dump(new)
    metadata = ['setup.py', 'setup.cfg', 'MANIFEST.in']
    assert all(baseline(p) == (ROOT/p).read_bytes() for p in metadata)
    wheel = args.dist_dir/'gem-0.1.12-py3-none-any.whl'
    sdist = args.dist_dir/'gem-0.1.12.tar.gz'
    with zipfile.ZipFile(wheel) as archive:
        assert all(archive.read(p) == (ROOT/p).read_bytes() for p in runtime)
    with tarfile.open(sdist) as archive:
        assert all(archive.extractfile('gem-0.1.12/'+p).read() == (ROOT/p).read_bytes() for p in runtime)
    smokes = {kind: json.loads((evidence/('phase4g4r-'+kind+'-smoke.json')).read_text())
              for kind in ('wheel', 'sdist')}
    assert all(s['status'] == 'passed' and s['original_regressions'] == 14 and s['isolated']
               for s in smokes.values())
    cross_modules = {'tests.test_geometry_audit', 'tests.test_geometry_repair',
                     'tests.test_angles_refraction', 'tests.test_planes', 'tests.test_numerical_robustness'}
    cross_cases = [c for c in full_cases if c['classname'] in cross_modules]
    assert len(cross_cases) == 580 and all(c['status'] == 'passed' for c in cross_cases)
    result = dict(base=BASE, status='passed', baseline=baseline_counts,
                  original_reproductions=original_counts, reproduction_cases=original_cases,
                  full_suite=full_counts, repaired_audit_cases=repaired,
                  new_independent_cases=160, cross_module_cases=len(cross_cases),
                  preservation=dict(original_audit_assertions_inputs_tolerances_unchanged=True,
                    strict_defect_markers_removed=removed, other_runtime_modules_unchanged=untouched,
                    public_signatures_unchanged=True, only_target_function_bodies_changed=True,
                    packaging_dependencies_metadata_unchanged=metadata),
                  installed_smokes=smokes, artifacts={p.name: dict(
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
                    runtime_modules_matching_source=len(runtime)) for p in (wheel, sdist)})
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print('Verified preservation, executed tests and both installed artifacts:', args.output)


if __name__ == '__main__':
    main()
