"""Reporting contracts use synthetic timings, never wall-clock assertions."""
import copy
import json
from pathlib import Path
import subprocess
import sys

import pytest

from benchmarks import compare as reporting


def artifact(latency=10, spread=(.99, 1, 1.01)):
    blocks = []
    for factor in spread:
        blocks.append({'seconds_per_operation': [latency*factor*1e-6]*7,
                       'operations_per_trial': 128, 'trials': 7})
    return {'schema': 'gem-benchmark-v1',
            'environment': dict.fromkeys(reporting.ENV_KEYS, 'same'),
            'methodology': dict.fromkeys(reporting.METHOD_KEYS, 'same'),
            'results': {'dot': {'definition_sha256': 'fixed-data-v1',
                                'unit': 'call', 'blocks': blocks}}}


@pytest.mark.parametrize('latency,status', [(20, 'possible_regression'),
    (5, 'possible_improvement'), (10.2, 'inconclusive'), (10, 'inconclusive')])
def test_independent_changes(latency, status):
    report = reporting.compare([artifact()], [artifact(latency)])
    row = report['results']['dot']
    assert row['status'] == status
    assert row['slowdown_ratio'] == pytest.approx(latency/10)
    assert row['speedup_ratio'] == pytest.approx(10/latency)
    assert row['delta_us'] == pytest.approx(latency-10)
    assert report['advisory_only'] is True


def test_large_variability_is_inconclusive():
    assert reporting.compare([artifact(spread=(.5, 1, 2))],
        [artifact(15, (.5, 1, 2))])['results']['dot']['status'] == 'inconclusive'


def test_within_trial_variability_is_reported():
    current = artifact(20)
    for block in current['results']['dot']['blocks']:
        block['seconds_per_operation'] = [1e-6, 2e-6, 3e-6, 20e-6, 40e-6, 50e-6, 60e-6]
    row = reporting.compare([artifact()], [current])['results']['dot']
    assert row['within_block_relative_mad']['current'] > .8
    assert row['status'] == 'inconclusive'


@pytest.mark.parametrize('key', reporting.ENV_KEYS)
def test_environment_mismatch_prevents_claim(key):
    current = artifact(20)
    current['environment'][key] = 'different'
    report = reporting.compare([artifact()], [current])
    assert report['results']['dot']['status'] == 'environment_not_comparable'
    assert 'environment.' + key in report['environment_mismatches']


def test_unknown_host_is_not_equivalence():
    baseline = artifact()
    del baseline['environment']['host_id']
    assert reporting.compare([baseline], [baseline])['results']['dot']['status'] == 'environment_not_comparable'


@pytest.mark.parametrize('key', reporting.METHOD_KEYS)
def test_methodology_mismatch(key):
    current = artifact(20)
    current['methodology'][key] = 'different'
    assert reporting.compare([artifact()], [current])['results']['dot']['status'] == 'environment_not_comparable'


@pytest.mark.parametrize('field,value', [('definition_sha256', 'new-data'),
    ('definition_version', 'v2'), ('metadata', {'count': 99}),
    ('unit', 'batch'), ('output_points', 100)])
def test_changed_workload(field, value):
    current = artifact()
    current['results']['dot'][field] = value
    assert reporting.compare([artifact()], [current])['results']['dot'] == {'status': 'workload_changed'}


def test_missing_and_new_workloads():
    current = artifact()
    current['results']['new'] = current['results'].pop('dot')
    results = reporting.compare([artifact()], [current])['results']
    assert results['dot']['status'] == 'missing_current'
    assert results['new']['status'] == 'missing_baseline'


@pytest.mark.parametrize('blocks,trials', [(2, 7), (3, 3)])
def test_insufficient_repetition(blocks, trials):
    current = artifact(20)
    current['results']['dot']['blocks'] = current['results']['dot']['blocks'][:blocks]
    for block in current['results']['dot']['blocks']:
        block['seconds_per_operation'] = block['seconds_per_operation'][:trials]
        block['trials'] = trials
    assert reporting.compare([artifact()], [current])['results']['dot']['status'] == 'insufficient_evidence'


@pytest.mark.parametrize('value', [float('nan'), float('inf'), 0, -1, True, 'fast'])
def test_invalid_samples(value):
    current = artifact()
    current['results']['dot']['blocks'][0]['seconds_per_operation'][0] = value
    with pytest.raises(ValueError):
        reporting.compare([artifact()], [current])


def test_trial_count_validation():
    current = artifact()
    current['results']['dot']['blocks'][0]['trials'] = 8
    with pytest.raises(ValueError):
        reporting.compare([artifact()], [current])


@pytest.mark.parametrize('value', [True, 0, -1, 1.5])
def test_operation_count_validation(value):
    current = artifact()
    current['results']['dot']['blocks'][0]['operations_per_trial'] = value
    with pytest.raises(ValueError):
        reporting.compare([artifact()], [current])


def test_absolute_noise_floor():
    row = reporting.compare([artifact(.01)], [artifact(.02)])['results']['dot']
    assert row['status'] == 'inconclusive'


def test_comparison_preserves_inputs_and_is_repeatable():
    before, after = artifact(), artifact(20)
    saved = copy.deepcopy((before, after))
    first = reporting.compare([before], [after])
    assert first == reporting.compare([before], [after])
    assert (before, after) == saved


def test_duplicate_artifact_is_not_repeated_evidence(tmp_path):
    path = tmp_path/'same.json'
    path.write_text(json.dumps(artifact()))
    data = reporting.load(path)
    with pytest.raises(ValueError, match='duplicate'):
        reporting.compare([data, data], [data])


def test_historical_provenance_is_confined_to_one_artifact(tmp_path):
    baseline, current = artifact(), artifact(5)
    historical = {'rounds': [{'dot': {'before': b, 'after': a}} for b, a in zip(
        baseline['results']['dot']['blocks'], current['results']['dot']['blocks'])]}
    path = tmp_path/'paired.json'
    path.write_text(json.dumps(historical))
    report = reporting.compare([reporting.load(path, 'before')], [reporting.load(path, 'after')])
    assert report['results']['dot']['status'] == 'possible_improvement'
    assert report['results']['dot']['resampling'] == 'paired blocks'
    other = tmp_path/'other.json'
    historical['extra'] = 'different artifact'
    other.write_text(json.dumps(historical))
    report = reporting.compare([reporting.load(path, 'before')], [reporting.load(other, 'after')])
    assert report['results']['dot']['status'] == 'workload_changed'


def test_cli_reports_regression_without_timing_gate(tmp_path):
    baseline, current, output = (tmp_path/name for name in ('before.json', 'after.json', 'report.json'))
    baseline.write_text(json.dumps(artifact()))
    current.write_text(json.dumps(artifact(20)))
    script = Path(__file__).resolve().parents[1]/'benchmarks'/'compare.py'
    result = subprocess.run([sys.executable, str(script), '--baseline', str(baseline),
        '--current', str(current), '--output', str(output)], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert 'possible_regression' in result.stdout
    assert json.loads(output.read_text())['counts'] == {'possible_regression': 1}


def test_workload_fingerprints_exclude_mathematical_implementation():
    from benchmarks import run_core
    calls, details = run_core.core.prepare()
    first = run_core.definitions(details)
    assert set(first) == set(calls)
    assert first == run_core.definitions(details)
    altered = copy.deepcopy(details)
    altered['vector3_dot']['domain'] = 'different inputs'
    changed = run_core.definitions(altered)
    assert first['vector3_dot'] != changed['vector3_dot']
    assert first['matrix4_inverse'] == changed['matrix4_inverse']


def test_quick_and_stress_partition_and_cli(tmp_path):
    from benchmarks import run_core
    calls, _ = run_core.core.prepare()
    assert run_core.STRESS <= calls.keys()
    script = Path(run_core.__file__)
    quick = subprocess.run([sys.executable, str(script), '--list'],
        capture_output=True, text=True, check=True).stdout.splitlines()
    extended = subprocess.run([sys.executable, str(script), '--suite', 'extended', '--list'],
        capture_output=True, text=True, check=True).stdout.splitlines()
    assert not set(quick) & set(extended)
    assert set(quick) | set(extended) == set(calls)
    output = tmp_path/'smoke.json'
    subprocess.run([sys.executable, str(script), '--case', 'vector3_dot', '--rounds', '1',
        '--trials', '3', '--target-seconds', '.0001', '--output', str(output)],
        capture_output=True, text=True, check=True)
    result = reporting.load(output)
    assert set(result['results']) == {'vector3_dot'}
    assert result['results']['vector3_dot']['blocks'][0]['trials'] == 3
    assert reporting.compare([result], [result])['results']['vector3_dot']['status'] == 'insufficient_evidence'


@pytest.mark.parametrize('options', [['--case', 'missing'],
    ['--case', 'vector3_dot', '--case', 'vector3_dot'], ['--trials', '0'],
    ['--rounds', '0'], ['--target-seconds', 'nan']])
def test_runner_rejects_invalid_measurement_options(options):
    script = Path(__file__).resolve().parents[1]/'benchmarks'/'run_core.py'
    result = subprocess.run([sys.executable, str(script), '--list'] + options,
                            capture_output=True, text=True)
    assert result.returncode == 2


def test_malformed_comparison_has_command_error(tmp_path):
    path = tmp_path/'invalid.json'
    path.write_text('{')
    script = Path(__file__).resolve().parents[1]/'benchmarks'/'compare.py'
    result = subprocess.run([sys.executable, str(script), '--baseline', str(path),
        '--current', str(path)], capture_output=True, text=True)
    assert result.returncode == 2


@pytest.mark.parametrize('relative,absolute', [(-1, .05), (.15, -1), (float('nan'), .05)])
def test_invalid_review_margins(relative, absolute):
    with pytest.raises(ValueError):
        reporting.compare([artifact()], [artifact()], relative, absolute)
