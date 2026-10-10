"""Independent negative checks for the documentation-only freeze and learning guard."""
import importlib.util
import hashlib
import json
from pathlib import Path
import subprocess

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def checker(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('docs_checker', ROOT/'tools/check_architecture_docs.py')
    tool = importlib.util.module_from_spec(spec); spec.loader.exec_module(tool)
    for name, text in {'gem/example.py': 'VALUE = 1\n', 'setup.py': '# metadata\n',
            'tests/existing.py': '# retained test\n',
            'docs/tutorials/vectors.md': '# Tutorial\n\n```python\nx = 1\nassert x == 1\n```\n'}.items():
        path = tmp_path/name; path.parent.mkdir(parents=True, exist_ok=True); path.write_bytes(text.encode())
    for args in (['init'], ['config', 'core.autocrlf', 'false'], ['add', '.'], ['-c', 'user.name=Fixture', '-c',
            'user.email=fixture@example.invalid', 'commit', '-m', 'fixture']):
        subprocess.run(['git', *args], cwd=tmp_path, capture_output=True, check=True)
    base = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=tmp_path, text=True).strip()
    monkeypatch.setattr(tool, 'ROOT', tmp_path); monkeypatch.setattr(tool, 'BASE', base)
    return tool, tmp_path


@pytest.mark.parametrize('name', ['gem/example.py', 'setup.py', 'tests/existing.py'])
def test_existing_scope_modifications_rejected(checker, name):
    tool, root = checker
    assert tool.check_scope()['protected_files_unchanged'] == 3
    (root/name).write_text('changed\n')
    with pytest.raises(ValueError, match='out-of-scope modification'):
        tool.check_scope()


def test_new_runtime_file_rejected(checker):
    tool, root = checker
    (root/'gem/extra.py').write_text('NEW_API = 1\n')
    with pytest.raises(ValueError, match='out-of-scope addition'):
        tool.check_scope()


def test_new_documentation_test_is_the_only_allowed_scoped_addition(checker):
    tool, root = checker
    (root/'tests/test_documentation_overhaul.py').write_text('# documentation checks\n')
    assert tool.check_scope()['new_documentation_test_modules'] == ['tests/test_documentation_overhaul.py']
    (root/'tests/unreviewed.py').write_text('# unexpected\n')
    with pytest.raises(ValueError, match='out-of-scope addition'):
        tool.check_scope()


def test_comment_only_example_edit_preserves_ast(checker):
    tool, root = checker; page = 'docs/tutorials/vectors.md'
    text = (root/page).read_text().replace('x = 1', '# Independent reference\nx = 1')
    (root/page).write_text(text)
    assert tool.preserve_learning_material([page]) == {
        'baseline_examples_preserved':1, 'tutorials_preserved':1, 'baseline_anchors_preserved':1}


@pytest.mark.parametrize('replacement', ['x = 2', ''])
def test_changed_or_removed_original_example_rejected(checker, replacement):
    tool, root = checker; page = 'docs/tutorials/vectors.md'
    (root/page).write_text('# Tutorial\n'+('```python\n'+replacement+'\n```' if replacement else ''))
    with pytest.raises(ValueError, match='removed/changed executable example'):
        tool.preserve_learning_material([page])


def test_old_heading_anchor_requires_compatible_destination(checker):
    tool, root = checker; page = 'docs/tutorials/vectors.md'
    text = (root/page).read_text().replace('# Tutorial', '# Renamed tutorial')
    (root/page).write_text(text)
    with pytest.raises(ValueError, match='removed historical anchor'):
        tool.preserve_learning_material([page])
    (root/page).write_text('<div id="tutorial"></div>\n\n'+text)
    assert tool.preserve_learning_material([page])['baseline_anchors_preserved'] == 1


@pytest.fixture
def wiki_update(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('wiki_update', ROOT/'tools/prepare_wiki.py')
    tool = importlib.util.module_from_spec(spec); spec.loader.exec_module(tool)
    wiki = tmp_path/'wiki'; wiki.mkdir()
    for name, text in {'Home.md':'Old portal\n', 'Historical-Home.md':'Historical text\n',
                       'Matrix-Class.md':'Historical mathematics\n'}.items():
        (wiki/name).write_text(text)
    for args in (['init'], ['config', 'core.autocrlf', 'false'], ['add', '.'], ['-c', 'user.name=Fixture', '-c',
            'user.email=fixture@example.invalid', 'commit', '-m', 'fixture']):
        subprocess.run(['git', *args], cwd=wiki, capture_output=True, check=True)
    package = tmp_path/'package'; pages = package/'pages'; pages.mkdir(parents=True)
    (pages/'Home.md').write_text('New navigation\n')
    manifest = {'observed_tip':subprocess.check_output(['git','rev-parse','HEAD'],cwd=wiki,text=True).strip(),
        'observed_files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in wiki.glob('*.md')},
        'files':{'Home.md':hashlib.sha256((pages/'Home.md').read_bytes()).hexdigest()},
        'preserve_existing':['Historical-Home.md','Matrix-Class.md']}
    (package/'manifest.json').write_text(json.dumps(manifest))
    monkeypatch.setattr(tool, 'PACKAGE', package)
    return tool, wiki, package, manifest


def test_reviewed_portal_update_preserves_history(wiki_update):
    tool, wiki, package, manifest = wiki_update
    assert tool.prepare(wiki) == 1 and (wiki/'Home.md').read_text() == 'Old portal\n'
    assert tool.prepare(wiki, True) == 1
    assert (wiki/'Home.md').read_text() == 'New navigation\n'
    for name in manifest['preserve_existing']:
        assert hashlib.sha256((wiki/name).read_bytes()).hexdigest() == manifest['observed_files'][name]


@pytest.mark.parametrize('field', ['observed_files','files'])
def test_portal_update_rejects_mismatches_before_writing(wiki_update, field):
    tool, wiki, package, manifest = wiki_update
    manifest[field]['Home.md'] = '0'*64
    (package/'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='differs|integrity mismatch'):
        tool.prepare(wiki, True)
    assert (wiki/'Home.md').read_text() == 'Old portal\n'
