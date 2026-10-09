"""Independent numerical and artifact checks for headless gallery examples."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import pytest
from examples.showcase.regenerate import generate, main
from examples.showcase.scenes import SCENES
from examples.showcase.verify import check_numerics, check_svg
from examples.showcase.svg import Canvas

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize('name', SCENES)
def test_diagram_structure_and_determinism(name):
    first, values = SCENES[name]()
    second, repeated = SCENES[name]()
    assert first == second and values == repeated
    assert check_svg(first) == [960, 560]
    assert len(ET.fromstring(first).findall('{http://www.w3.org/2000/svg}text')) >= 8


def test_independent_geometry_references():
    data = {name: scene()[1] for name, scene in SCENES.items()}
    result = check_numerics(data)
    assert result['cubic_dense_reference_max_distances'][1] < .05
    data['transforms']['results'][4][0][0] += 1
    with pytest.raises(AssertionError):
        check_numerics(data)


def test_committed_artifacts_match_generation():
    for name, payload in generate().items():
        assert (ROOT/'examples/showcase/output'/name).read_bytes() == payload


def test_output_scope_repeated_runs_and_verification(tmp_path):
    sentinel = tmp_path/'unrelated.txt'
    sentinel.write_text('preserve')
    main(['--output-dir', str(tmp_path)])
    initial = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    main(['--output-dir', str(tmp_path)])
    main(['--output-dir', str(tmp_path), '--verify'])
    assert initial == {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    (tmp_path/'vectors.svg').write_text('altered')
    with pytest.raises(SystemExit, match='artifact mismatch'):
        main(['--output-dir', str(tmp_path), '--verify'])
    assert sentinel.read_text() == 'preserve'


def test_manifest_and_svg_escaping():
    report = json.loads(generate()['measurements.json'])
    assert report['randomness'] == 'none'
    assert report['scenes']['bezier']['quadratic_midpoint'] == [2., 2.]
    c = Canvas('A < B & C', 'quoted "description"')
    c.text(10, 10, '<script> & text')
    parsed = ET.fromstring(c.finish())
    assert parsed.find('{http://www.w3.org/2000/svg}title').text == 'A < B & C'
    assert '&lt;script&gt;' in c.finish()
