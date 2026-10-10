"""Release metadata and narrow packaging-freeze guards."""
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import sys

import pytest

ROOT=Path(__file__).resolve().parents[1]


def test_authoritative_version_and_dependency_free_package_import():
    from gem import __version__
    assert __version__ == '1.0.0'
    tree=ast.parse((ROOT/'gem/_version.py').read_text())
    assert not any(isinstance(n,(ast.Import,ast.ImportFrom)) for n in ast.walk(tree))
    config=(ROOT/'pyproject.toml').read_text()
    assert 'dynamic = ["version"]' in config
    assert 'version = {attr = "gem._version.__version__"}' in config
    assert not re.search(r'^version\s*=\s*"', config, re.M)


def test_package_identity_and_explicit_runtime_dependency():
    config=(ROOT/'pyproject.toml').read_text()
    assert 'name = "gem"' in config
    assert 'dependencies = ["six"]' in config
    assert 'requires-python = ">=3.10"' in config
    assert 'packages = ["gem", "gem.experimental"]' in config
    assert 'Python :: 2' not in config
    assert 'BSD-2-Clause' in config
    assert 'explosiveduck' not in config


def test_metadata_changes_do_not_modify_frozen_mathematics():
    import subprocess
    base='15fbce6'
    for path in (ROOT/'gem').rglob('*.py'):
        name=str(path.relative_to(ROOT))
        if name in ('gem/__init__.py','gem/_version.py'): continue
        assert path.read_bytes() == subprocess.check_output(['git','show',base+':'+name],cwd=ROOT)


@pytest.mark.parametrize('name',['gem/__init__.py','gem/_version.py','setup.py','pyproject.toml'])
def test_packaging_manifest_rejects_unreviewed_content(tmp_path,monkeypatch,name):
    spec=importlib.util.spec_from_file_location('scope_guard',ROOT/'tools/check_architecture_docs.py')
    tool=importlib.util.module_from_spec(spec);spec.loader.exec_module(tool)
    import subprocess
    (tmp_path/'gem').mkdir();(tmp_path/'gem/original.py').write_text('VALUE=1\n')
    subprocess.run(['git','init'],cwd=tmp_path,check=True,capture_output=True)
    subprocess.run(['git','add','.'],cwd=tmp_path,check=True)
    subprocess.run(['git','-c','user.name=Fixture','-c','user.email=f@example.invalid','commit','-m','fixture'],cwd=tmp_path,check=True,capture_output=True)
    base=subprocess.check_output(['git','rev-parse','HEAD'],cwd=tmp_path,text=True).strip()
    path=tmp_path/name;path.parent.mkdir(exist_ok=True);path.write_text('reviewed\n')
    (tmp_path/'tools').mkdir()
    (tmp_path/'tools/phase5a-packaging-scope.json').write_text(json.dumps({'base':base,'files':{name:hashlib.sha256(path.read_bytes()).hexdigest()}}))
    monkeypatch.setattr(tool,'ROOT',tmp_path);monkeypatch.setattr(tool,'BASE',base)
    assert tool.check_scope()['reviewed_packaging_files']==[name]
    path.write_text('unreviewed\n')
    with pytest.raises(ValueError,match='unreviewed packaging change'): tool.check_scope()
