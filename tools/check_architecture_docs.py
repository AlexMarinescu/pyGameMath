"""Check architecture links, source declarations, examples and unchanged scope.

Standard-library documentation tooling, not an installed gem runtime module.
Run from any working directory: python tools/check_architecture_docs.py --examples.
"""
import argparse
import ast
import hashlib
import importlib
import inspect
import json
from pathlib import Path
import re
import subprocess
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
BASE = '15fbce64fa87008203890142070dbf30ea802ebc'
PAGES = ['ROADMAP.md', 'docs/README.md', 'docs/development/roadmap.md',
         'docs/development/verification.md'] + [
    'docs/architecture/' + name + '.md' for name in
    ('philosophy', 'conventions', 'compatibility', 'api-inventory')]


def declarations():
    result = {}
    for source in sorted((ROOT/'gem').glob('*.py')):
        rows = []
        for node in ast.parse(source.read_text()).body:
            if isinstance(node, ast.FunctionDef) and not node.name.startswith('_'):
                rows.append(node.name+'('+ast.unparse(node.args)+')')
            elif isinstance(node, ast.ClassDef) and not node.name.startswith('_'):
                rows.append(node.name)
                for member in node.body:
                    if isinstance(member, ast.FunctionDef) and (
                        not member.name.startswith('_') or member.name.startswith('__')):
                        rows.append(node.name+'.'+member.name+'('+ast.unparse(member.args)+')')
                    elif isinstance(member, ast.Assign) and isinstance(member.value, ast.Name):
                        for target in member.targets:
                            if isinstance(target, ast.Name) and target.id.startswith('__'):
                                rows.append(node.name+'.'+target.id+' = '+node.name+'.'+member.value.id)
        if rows:
            result['.'.join(source.relative_to(ROOT).with_suffix('').parts)] = rows
    return result


def anchors(text):
    result = set(re.findall(r'\bid=["\']([^"\']+)["\']', text))
    seen = {}
    for title in re.findall(r'^#{1,6}\s+(.+?)\s*#*$', text, re.MULTILINE):
        slug = re.sub(r'[^\w\- ]', '', title.replace('`', '').lower()).replace(' ', '-')
        occurrence = seen.get(slug, 0)
        seen[slug] = occurrence+1
        result.add(slug + ('-'+str(occurrence) if occurrence else ''))
    return result


def check_scope():
    """Freeze runtime, existing tests and packaging at the reviewed master.

    Baseline changes require review; HEAD is deliberately not the reference.
    Reviewed packaging changes are admitted only by exact content fingerprints.
    Runtime mathematics remains frozen; arbitrary additions remain rejected.
    """
    paths = subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE],
                                    cwd=ROOT, text=True).splitlines()
    def protected(name):
        return name.startswith(('gem/', 'tests/', 'examples/', 'benchmarks/')) or name in (
            'setup.py', 'setup.cfg', 'MANIFEST.in', 'LICENSE', 'README.rst',
            'requirements-audit.txt', 'requirements-docs.txt', 'requirements-docs-browser.txt',
            '.travis.yml')
    originals = {name for name in paths if protected(name)}
    current = subprocess.check_output(['git', 'ls-files', '--cached', '--others',
        '--exclude-standard'], cwd=ROOT, text=True).splitlines()
    additions = {name for name in current if protected(name)}-originals
    ci_reviewed = {}
    ci_manifest = ROOT/'tools/phase5b-ci-scope.json'
    if ci_manifest.exists():
        ci_scope = json.loads(ci_manifest.read_text())
        if ci_scope['base'] != BASE:
            raise ValueError('CI scope baseline mismatch')
        ci_reviewed = ci_scope['files']
        if set(ci_reviewed)-{'tests/test_documentation_overhaul.py',
                             'tests/test_release_packaging.py', 'benchmarks/run_core.py'}:
            raise ValueError('mathematics cannot be exempted by CI scope')
        for name, digest in ci_reviewed.items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
                raise ValueError('unreviewed CI portability change: '+name)
    reviewed = {}
    manifest = ROOT/'tools/phase5a-packaging-scope.json'
    if manifest.exists():
        scope = json.loads(manifest.read_text())
        if scope['base'] != BASE:
            raise ValueError('packaging scope baseline mismatch')
        reviewed = scope['files']
        permitted = {'pyproject.toml', 'setup.py', 'setup.cfg', 'MANIFEST.in',
                     'README.rst', 'requirements-build.txt', 'gem/__init__.py',
                     'gem/_version.py', 'tests/test_core_packaging.py',
                     'tests/test_release_packaging.py'}
        if set(reviewed)-permitted:
            raise ValueError('mathematics cannot be exempted by packaging scope')
        for name, digest in reviewed.items():
            if name not in ci_reviewed and hashlib.sha256((ROOT/name).read_bytes()).hexdigest() != digest:
                raise ValueError('unreviewed packaging change: '+name)
    allowed = {'tests/test_documentation_overhaul.py',
               'tests/test_release_workflows.py'} | set(reviewed)
    if additions-allowed:
        raise ValueError('out-of-scope addition: '+', '.join(sorted(additions-allowed)))
    fingerprint = hashlib.sha256()
    for name in sorted(originals):
        content = (ROOT/name).read_bytes()
        if name not in reviewed and name not in ci_reviewed and content != subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT):
            raise ValueError('out-of-scope modification: '+name)
        fingerprint.update(name.encode()); fingerprint.update(content)
    return {'protected_files_unchanged': len(originals-set(reviewed)-set(ci_reviewed)),
            'protected_sha256': fingerprint.hexdigest(),
            'new_documentation_test_modules': sorted(additions),
            'reviewed_ci_files': sorted(ci_reviewed),
            'reviewed_packaging_files': sorted(reviewed)}


def preserve_learning_material(pages):
    """Keep every baseline executable block in order, permitting comments only."""
    baseline = set(subprocess.check_output(['git', 'ls-tree', '-r', '--name-only', BASE],
                                         cwd=ROOT, text=True).splitlines())
    blocks = tutorials = old_anchors = 0
    for name in pages:
        if name not in baseline:
            continue
        previous = subprocess.check_output(['git', 'show', BASE+':'+name], cwd=ROOT, text=True)
        current = (ROOT/name).read_text()
        required = anchors(re.sub(r'```.*?```', '', previous, flags=re.DOTALL))
        present = anchors(re.sub(r'```.*?```', '', current, flags=re.DOTALL))
        if required-present:
            raise ValueError('removed historical anchor: '+name+': '+', '.join(sorted(required-present)))
        old_anchors += len(required)
        old = re.findall(r'```python\n(.*?)```', previous, re.DOTALL)
        new = re.findall(r'```python\n(.*?)```', current, re.DOTALL)
        expected = [ast.dump(ast.parse(code)) for code in old]
        actual = [ast.dump(ast.parse(code)) for code in new]
        cursor = 0
        for block in expected:
            while cursor < len(actual) and actual[cursor] != block:
                cursor += 1
            if cursor == len(actual):
                raise ValueError('removed/changed executable example: '+name)
            cursor += 1
            blocks += 1
        if name.startswith('docs/tutorials/') and old:
            tutorials += 1
    return {'baseline_examples_preserved': blocks, 'tutorials_preserved': tutorials,
            'baseline_anchors_preserved': old_anchors}


def check_api_reference(declared, runtime=False):
    """Check declaration-to-reference coverage, constants and installed signatures."""
    coverage = json.loads((ROOT/'docs/api/coverage.json').read_text())
    entries = coverage['declarations']
    if set(entries) != set(declared) or coverage['omitted_declarations']:
        raise ValueError('API reference module coverage mismatch')
    count = 0
    for module, rows in declared.items():
        if [entry['signature'] for entry in entries[module]] != rows:
            raise ValueError('API reference declaration coverage mismatch: '+module)
        imported = importlib.import_module(module) if runtime else None
        for entry in entries[module]:
            signature = entry['signature']
            page = (ROOT/entry['page']).resolve()
            if not page.is_relative_to(ROOT/'docs/api'):
                raise ValueError('API coverage page escapes reference directory')
            text = page.read_text()
            if text.count('| `'+signature+'` |') != 1:
                raise ValueError('missing/duplicate API declaration row: '+signature)
            description = text.split('| `'+signature+'` |', 1)[1].split('\n', 1)[0].strip(' |')
            if not description:
                raise ValueError('empty API behavior description: '+signature)
            if runtime:
                path = signature.split('(')[0].split(' = ')[0]
                value = imported
                for part in path.split('.'):
                    value = getattr(value, part)
                if ' = ' in signature:
                    target = imported
                    for part in signature.split(' = ')[1].split('.'):
                        target = getattr(target, part)
                    if value is not target:
                        raise ValueError('runtime alias mismatch: '+signature)
                elif '(' in signature:
                    expected = signature[signature.index('(')+1:-1]
                    actual = ast.unparse(ast.parse('def check'+str(inspect.signature(value))+': pass').body[0].args)
                    if actual != expected:
                        raise ValueError('runtime signature mismatch: '+module+'.'+signature)
                elif not inspect.isclass(value):
                    raise ValueError('runtime class mismatch: '+signature)
            count += 1
    discovered = {}
    for source in sorted((ROOT/'gem').glob('*.py')):
        names = []
        for node in ast.parse(source.read_text()).body:
            if isinstance(node, ast.Assign):
                names.extend(target.id for target in node.targets
                             if isinstance(target, ast.Name) and not target.id.startswith('_'))
        if names:
            discovered['gem.'+source.stem] = names
    if discovered != coverage['constants']:
        raise ValueError('exposed API constant coverage mismatch')
    if runtime:
        import ctypes
        common = importlib.import_module('gem.common')
        if common.GLfloat is not ctypes.c_float:
            raise ValueError('runtime GLfloat mismatch')
        vector = importlib.import_module('gem.vector')
        for name in coverage['constants']['gem.vector']:
            size = int(name[-1])
            expected = [1.0 if name.startswith('I') else 0.0]*size
            if getattr(vector, name) != expected:
                raise ValueError('runtime reference buffer mismatch: '+name)
        for source, targets in coverage['reexports'].items():
            if isinstance(targets, str):
                targets = [targets]
                shim = importlib.import_module(source.rsplit('.', 1)[0])
                names = [source.rsplit('.', 1)[1]]
            else:
                shim = importlib.import_module(source)
                names = [target.rsplit('.', 1)[1] for target in targets]
            for name, target in zip(names, targets):
                module, symbol = target.rsplit('.', 1)
                if getattr(shim, name) is not getattr(importlib.import_module(module), symbol):
                    raise ValueError('runtime compatibility identity mismatch: '+source)
    return {'declarations': count, 'constants': sum(map(len, discovered.values())),
            'runtime_signatures_checked': runtime, 'omitted_declarations': []}


def check(examples=False, getting_started=False, package_root=None, api_reference=False, tutorials=False, showcase=False, website=False):
    pages = PAGES + (['README.md'] + ['docs/getting-started/'+name+'.md'
                     for name in ('README', 'installation', 'quick-start', 'verification')]
                     if getting_started else [])
    if api_reference:
        pages += [str(page.relative_to(ROOT)) for page in sorted((ROOT/'docs/api').glob('*.md'))]
    if tutorials:
        pages += [str(page.relative_to(ROOT)) for page in sorted((ROOT/'docs/tutorials').glob('*.md'))]
    if showcase:
        pages += [str(page.relative_to(ROOT)) for page in sorted((ROOT/'docs/examples').rglob('*.md'))]
        pages += ['examples/showcase/README.md']
    if website:
        pages += ['docs/index.md', 'docs/architecture/notation.md', 'docs/wiki-migration/README.md']
        pages += ['docs/development/'+name+'.md' for name in ('documentation-stack', 'contributing', 'releases', 'website', 'website-verification', 'wiki', 'visual-overhaul')]
    report = {'base': BASE, 'pages': pages, 'links_checked': 0,
              'source_declarations_checked': 0, 'examples_executed': 0, 'examples_by_page': {}}
    inventory = (ROOT/'docs/architecture/api-inventory.md').read_text()
    declared = declarations()
    for module, rows in declared.items():
        heading = '### '+module+' declarations'
        if heading not in inventory:
            raise ValueError('missing module catalog: '+module)
        section = inventory.split(heading, 1)[1].split('\n### ', 1)[0]
        documented = re.findall(r'^\| (?:Function|Class|Method|Alias) \| `([^`]+)` \|$',
                                section, re.MULTILINE)
        if documented != rows:
            raise ValueError('source/catalog mismatch: '+module)
        report['source_declarations_checked'] += len(rows)
    if package_root is None:
        sys.path.insert(0, str(ROOT))
    else:
        package_root = package_root.resolve()
        sys.path.insert(0, str(package_root))
        import gem
        if not Path(gem.__file__).resolve().is_relative_to(package_root):
            raise ValueError('examples imported gem outside installed package root')
        report['example_package'] = str(Path(gem.__file__).resolve())
    if api_reference:
        report['api_reference'] = check_api_reference(declared, runtime=examples)
    report.update(preserve_learning_material(pages))
    for relative in pages:
        page = ROOT/relative
        text = page.read_text()
        # Exclude fenced code before interpreting Markdown links/headings.
        prose = re.sub(r'```.*?```', '', text, flags=re.DOTALL)
        for target in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)', prose):
            parts = urlsplit(target)
            if parts.scheme or parts.netloc:
                continue
            dest = (page.parent/unquote(parts.path)).resolve() if parts.path else page
            if not dest.is_relative_to(ROOT) or not dest.exists():
                raise ValueError(relative+': missing/internal escape link '+target)
            if parts.fragment:
                if not dest.is_file() or unquote(parts.fragment) not in anchors(dest.read_text()):
                    raise ValueError(relative+': missing anchor '+target)
            report['links_checked'] += 1
        if examples:
            for index, code in enumerate(re.findall(r'```python\n(.*?)```', text, re.DOTALL)):
                exec(compile(code, relative+':example'+str(index+1), 'exec'), {})
                report['examples_executed'] += 1
                report['examples_by_page'][relative] = report['examples_by_page'].get(relative, 0) + 1
    report.update(check_scope())
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--examples', action='store_true')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--getting-started', action='store_true')
    parser.add_argument('--api-reference', action='store_true')
    parser.add_argument('--tutorials', action='store_true')
    parser.add_argument('--showcase', action='store_true')
    parser.add_argument('--website', action='store_true')
    parser.add_argument('--package-root', type=Path,
                        help='execute examples against this installed site-packages directory')
    args = parser.parse_args()
    try:
        report = check(args.examples, args.getting_started, args.package_root, args.api_reference, args.tutorials, args.showcase, args.website)
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        parser.error(str(error))
    rendered = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered)
    print(rendered)


if __name__ == '__main__':
    main()
