"""Documentation safeguards use standard-library tools, not runtime dependencies."""
import importlib.util
import json
from pathlib import Path
import subprocess
import pytest

ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT/'tools'/(name+'.py'))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def site_fixture(root):
    required = ['api/vector', 'api/matrix', 'api/quaternion', 'api/spherical-harmonics',
                'examples', 'development/canonical-roadmap', 'architecture/notation']
    for name in required:
        p = root/name/'index.html'; p.parent.mkdir(parents=True,exist_ok=True)
        p.write_text('<html><body id="top">'+('<math></math>'*16+'<span class="kn">import</span>' if name.endswith('notation') else '')+'</body></html>')
    root.joinpath('index.html').write_text('<html>'+''.join('<a class="md-nav__link" href="'+p+'/">entry</a>' for p in required)+'</html>')
    search = root/'search/search_index.json';search.parent.mkdir()
    search.write_text(json.dumps({'docs':[{'location':'api/vector/','text':'rotate_coefficients quat_slerp findDrawingPoints viewport'}]}))
    return root


def test_built_site_validation_and_broken_fragment(tmp_path):
    tool = module('check_site');site_fixture(tmp_path)
    assert tool.check(tmp_path)['mathml_equations'] == 16
    p = tmp_path/'index.html';p.write_text(p.read_text()+'<a href="api/vector/#missing">bad</a>')
    with pytest.raises(ValueError,match='missing fragment'):
        tool.check(tmp_path)


def test_missing_asset_and_remote_css_rejected(tmp_path):
    tool = module('check_site');site_fixture(tmp_path)
    p = tmp_path/'index.html';p.write_text(p.read_text()+'<img src="absent.svg" alt="diagram">')
    with pytest.raises(ValueError,match='missing link'):
        tool.check(tmp_path)
    p.write_text(p.read_text().replace('<img src="absent.svg" alt="diagram">',''))
    (tmp_path/'style.css').write_text('a{background:url(https://example.invalid/image.png)}')
    with pytest.raises(ValueError,match='non-local CSS'):
        tool.check(tmp_path)


def test_wiki_guard_preserves_history_and_never_partially_writes(tmp_path, monkeypatch):
    tool = module('prepare_wiki')
    wiki=tmp_path/'wiki';wiki.mkdir();(wiki/'Home.md').write_text('Historical home\n');(wiki/'Old.md').write_text('Historical class\n')
    for args in [['init'],['add','.'],['-c','user.name=Fixture','-c','user.email=fixture@example.invalid','commit','-m','fixture']]:
        subprocess.run(['git','-C',str(wiki),*args],check=True,capture_output=True)
    tip=subprocess.check_output(['git','-C',str(wiki),'rev-parse','HEAD'],text=True).strip()
    package=tmp_path/'package';pages=package/'pages';pages.mkdir(parents=True)
    (pages/'Home.md').write_text('Canonical navigation\n');(pages/'Historical-Home.md').write_text('Historical home\n')
    import hashlib
    manifest={'observed_tip':tip,'preserve_existing':['Old.md'],'files':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in pages.iterdir()}}
    (package/'manifest.json').write_text(json.dumps(manifest));monkeypatch.setattr(tool,'PACKAGE',package)
    assert tool.prepare(wiki) == 2 and (wiki/'Home.md').read_text()=='Historical home\n'
    (pages/'Historical-Home.md').write_text('tampered')
    with pytest.raises(ValueError,match='historical Home differs'):
        tool.prepare(wiki,True)
    assert (wiki/'Home.md').read_text()=='Historical home\n'
    (pages/'Historical-Home.md').write_text('Historical home\n')
    tool.prepare(wiki,True)
    assert (wiki/'Old.md').read_text()=='Historical class\n'
    assert (wiki/'Historical-Home.md').read_text()=='Historical home\n'
    with pytest.raises(ValueError,match='clean'):
        tool.prepare(wiki)
