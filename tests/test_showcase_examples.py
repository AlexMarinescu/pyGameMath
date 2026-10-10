"""Independent numerical and artifact checks for headless gallery examples."""
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import pytest
from examples.showcase.regenerate import generate, main
from examples.showcase.scenes import SCENES
from examples.showcase.verify import check_numerics, check_svg, check_artifact
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
        check_artifact(name, (ROOT/'examples/showcase/output'/name).read_bytes(), payload)


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


def encode_measurements(report):
    return (json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()


def test_audited_libm_variation_preserves_artifact_verification(monkeypatch):
    import math
    # Replay the exact first differing libm operation, not a fabricated fixture.
    original = math.sin
    def windows_sine(angle):
        if angle.hex() == '0x1.921fb54442d18p-1':
            return float.fromhex('0x1.6a09e667f3bcdp-1')
        return original(angle)
    monkeypatch.setattr(math, 'sin', windows_sine)
    for name, payload in generate().items():
        check_artifact(name, (ROOT/'examples/showcase/output'/name).read_bytes(), payload)


@pytest.mark.parametrize('mutation', ['excess', 'unlisted', 'metadata', 'type', 'key', 'nan', 'serialization'])
def test_artifact_portability_policy_rejects_corruption(mutation):
    import math
    reference = (ROOT/'examples/showcase/output/measurements.json').read_bytes()
    report = json.loads(reference)
    frames = report['scenes']['quaternions']['frames']
    if mutation == 'excess': frames[1]['axes'][0][0] += 2*math.ulp(1.)
    elif mutation == 'unlisted': frames[2]['quaternion'][0] = math.nextafter(frames[2]['quaternion'][0], 1.)
    elif mutation == 'metadata': report['randomness'] = 'changed'
    elif mutation == 'type': frames[1]['axes'][0][0] = str(frames[1]['axes'][0][0])
    elif mutation == 'key': del frames[1]['axes']
    payload = encode_measurements(report)
    if mutation == 'nan': payload = payload.replace(b'0.8660254037844387', b'NaN', 1)
    elif mutation == 'serialization': payload = payload.rstrip()
    with pytest.raises((AssertionError, KeyError, TypeError)):
        check_artifact('measurements.json', reference, payload)


def test_portability_allowance_cannot_hide_shared_mathematical_error():
    reference = json.loads((ROOT/'examples/showcase/output/measurements.json').read_bytes())
    reference['scenes']['quaternions']['frames'][1]['axes'][0][0] += 1e-14
    payload = encode_measurements(reference)
    with pytest.raises(AssertionError, match='closed-form'):
        check_artifact('measurements.json', payload, payload)


def test_svg_verification_remains_byte_exact():
    reference = generate()['quaternions.svg']
    with pytest.raises(AssertionError):
        check_artifact('quaternions.svg', reference, reference.replace(b'120', b'121'))
