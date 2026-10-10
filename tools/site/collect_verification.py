"""Collect executed Phase 4F-C results and assert documentation-only scope."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[2]
BASE = '3e714fe949b7a6b7724d5c0da3395ee92483265f'


def sha(data):
    return hashlib.sha256(data).hexdigest()


def junit(name):
    path = Path('/tmp/phase4fc-'+name+'.xml')
    root = ET.parse(path).getroot()
    groups = {key: [] for key in ('passed','failed','xfail','skipped','errors')}
    for case in root.findall('.//testcase'):
        key = 'passed'
        if case.find('failure') is not None: key = 'failed'
        elif case.find('error') is not None: key = 'errors'
        elif case.find('skipped') is not None:
            key = 'xfail' if case.find('skipped').get('type') == 'pytest.xfail' else 'skipped'
        groups[key].append(case.get('classname')+'::'+case.get('name'))
    return {**{key: len(value) for key,value in groups.items()},
        'total':sum(map(len,groups.values())), 'junit_sha256':sha(path.read_bytes()),
        'elapsed_seconds':sum(float(s.get('time')) for s in root.findall('testsuite')),
        'new_documentation_tests':[n for n in groups['passed'] if n.startswith('tests.test_documentation_overhaul::')]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    tests = {kind:junit(kind) for kind in ('focused','full')}
    assert tests['full']['passed'] == tests['full']['total'] == 3959
    assert tests['focused']['passed'] == tests['focused']['total'] == 15
    assert len(tests['full']['new_documentation_tests']) == 12
    docs = {}
    for kind in ('source','wheel','sdist'):
        suffix = 'docs' if kind == 'source' else kind+'-docs'
        data = json.loads(Path('/tmp/phase4fc-'+suffix+'.json').read_text())
        assert data['source_declarations_checked'] == data['api_reference']['declarations'] == 268
        assert data['examples_executed'] == 43 and data['baseline_examples_preserved'] == 42
        assert data['tutorials_preserved'] == 9 and data['baseline_anchors_preserved'] == 219
        assert data['api_reference']['runtime_signatures_checked']
        docs[kind] = data
    site = json.loads(Path('/tmp/phase4fc-site.json').read_text())
    assert site['html_pages'] == 66 and site['external_asset_dependencies'] == 0
    browser = json.loads(Path('/tmp/phase4fc-browser/results.json').read_text())
    assert len(browser['layouts_checked']) == 56 and browser['mathml_notation_equations'] == 24
    assert not any(browser[k] for k in ('external_requests','javascript_errors','missing_resources'))
    assert min(browser['link_contrast_ratios'].values()) >= 4.5
    assert min(browser['brand_contrast_ratios'].values()) >= 4.5
    tracked = subprocess.check_output(['git','ls-tree','-r','--name-only',BASE],cwd=ROOT,text=True).splitlines()
    changed = []
    for name in tracked:
        old = subprocess.check_output(['git','show',BASE+':'+name],cwd=ROOT)
        if (ROOT/name).read_bytes() != old:
            assert name.startswith(('docs/','tools/')) or name in ('mkdocs.yml','.github/workflows/publish-docs.yml'), name
            changed.append(name)
    runtime = {str(p.relative_to(ROOT)):sha(p.read_bytes()) for p in (ROOT/'gem').rglob('*.py')}
    assert len(runtime) == 17
    installed = {}
    for kind in ('wheel','sdist'):
        package_root = Path(docs[kind]['example_package']).parent
        for name,digest in runtime.items():
            assert sha((package_root/Path(name).relative_to('gem')).read_bytes()) == digest
        installed[kind] = {'runtime_files_byte_identical':17,
            'package_root':str(package_root),
            'environment':'clean Phase 4G-5R installation; runtime unchanged by merged PR #59 and this phase'}
    wiki = Path('/tmp/phase4fc-wiki')
    manifest = json.loads((ROOT/'docs/wiki-migration/manifest.json').read_text())
    modified = subprocess.check_output(['git','-C',str(wiki),'diff','--name-only'],text=True).splitlines()
    assert sorted(modified) == sorted(manifest['files'])
    for name in manifest['preserve_existing']:
        assert sha((wiki/name).read_bytes()) == manifest['observed_files'][name]
    for name in manifest['files']:
        assert (wiki/name).read_bytes() == (ROOT/'docs/wiki-migration/pages'/name).read_bytes()
    plot = ROOT/'docs/assets/diagrams/legendre.svg'
    original = plot.read_bytes()
    subprocess.run([sys.executable,str(ROOT/'tools/site/legendre_diagram.py')],check=True,capture_output=True)
    assert plot.read_bytes() == original
    # Keep the runtime-only package metadata and every original test unchanged.
    assert not any(name.startswith(('gem/','tests/','examples/','benchmarks/')) for name in changed)
    screenshots = {p.name:sha(p.read_bytes()) for p in (ROOT/'audit/phase4fc-screenshots').glob('*.png')}
    assert set(screenshots) == {'home-light-viewport.png','home-dark-viewport.png','mobile-drawer.png'}
    for name, digest in screenshots.items():
        assert sha((Path('/tmp/phase4fc-browser')/name).read_bytes()) == digest
    result = {'schema':1,'phase':'4F-C','base':BASE,'branch':'docs/phase4fc-visual-overhaul',
        'environment':{'python':sys.version,'platform':platform.platform(),'other_python_versions_verified':False},
        'tests':tests,'documentation':docs,'built_site':site,'browser':browser,'installed':installed,
        'preservation':{'changed_original_files':changed,'all_other_original_files_unchanged':len(tracked)-len(changed),
            'runtime_sha256':runtime,'production_mathematics_unchanged':True,
            'public_api_signatures_unchanged':True,'metadata_versions_dependencies_unchanged':True,
            'existing_tests_examples_benchmarks_unchanged':True,'baseline_example_calculations_preserved':42,
            'baseline_tutorials_preserved':9,'baseline_heading_anchors_preserved':219},
        'legendre_plot':{'sha256':sha(original),'regeneration_byte_identical':True,
            'evaluates_core_api':True,'independent_analytical_checks':3},
        'wiki':{'observed_tip':manifest['observed_tip'],'existing_portal_published':True,
            'prepared_update_published':False,'reviewed_changed_files':modified,
            'other_pages_preserved':len(manifest['preserve_existing']),
            'website_url_source':manifest['website_url_source'],
            'live_website_check':'HTTP proxy rejected GitHub Pages domain with CONNECT 403; no bypass attempted'},
        'screenshots':screenshots,
        'commands_report':'audit/PHASE4FC-DOCUMENTATION.md'}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Collected tests, API/examples, 56 browser layouts, scope and Wiki preservation evidence')


if __name__ == '__main__':
    main()
