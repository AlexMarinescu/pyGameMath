"""Versioned current-core benchmarks; measurement tools are not gem runtime."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time

if __package__:
    from . import core_baseline as core
else:
    import core_baseline as core

STRESS = {'bezier_subdivide_1e-08', 'bezier_subdivide_1e-15',
          'sh_project_n1024_b1', 'sh_project_n1024_b3', 'sh_project_n1024_b5'}
SCHEMA = 'gem-benchmark-v1'
METHOD = {'timer': 'perf_counter', 'gc': 'timeit disables GC',
          'allocation': 'returned results included', 'execution': 'warm calibrated calls'}


def environment():
    cpu = 'unavailable'
    if Path('/proc/cpuinfo').exists():
        cpu = next((line.split(':', 1)[1].strip() for line in
                    Path('/proc/cpuinfo').read_text().splitlines()
                    if line.startswith('model name')), cpu)
    machine = Path('/etc/machine-id')
    host = (machine.read_text().strip() if machine.exists() else '') + platform.node()
    return {'python': sys.version, 'implementation': platform.python_implementation(),
            'platform': platform.platform(), 'architecture': platform.machine(),
            'cpu': cpu, 'cpu_count': os.cpu_count(),
            'host_id': hashlib.sha256(host.encode()).hexdigest() if host else None,
            'affinity': sorted(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else None,
            'six': importlib.metadata.version('six')}


def definitions(details):
    source = Path(core.__file__).read_text()
    dataset = source[source.index('def V('):source.index('def measure(')]
    return {name: hashlib.sha256(json.dumps(
        [dataset, 'core-v1', name, metadata, 'named workload call'],
        sort_keys=True).encode()).hexdigest() for name, metadata in details.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path)
    parser.add_argument('--suite', choices=('quick', 'extended'), default='quick')
    parser.add_argument('--case', action='append', default=[])
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--trials', type=int, default=7)
    parser.add_argument('--target-seconds', type=float, default=.01)
    parser.add_argument('--profiles', action='store_true')
    args = parser.parse_args()
    if args.rounds < 1 or args.trials < 3 or not math.isfinite(args.target_seconds) or args.target_seconds <= 0:
        parser.error('require positive rounds, >=3 trials and finite positive target')
    start = time.perf_counter()
    calls, details = core.prepare()
    setup = time.perf_counter() - start
    if set(args.case) - calls.keys():
        parser.error('unknown workload: ' + ', '.join(sorted(set(args.case) - calls.keys())))
    if len(args.case) != len(set(args.case)):
        parser.error('duplicate --case would repeat the same measured block')
    names = args.case or [name for name in calls if (name in STRESS) == (args.suite == 'extended')]
    if args.list:
        print('\n'.join(names))
        return
    if args.output is None:
        parser.error('--output is required unless --list is used')
    hashes = definitions(details)
    results = {name: {'definition_sha256': hashes[name], 'definition_version': 'core-v1',
                     'metadata': details[name], 'unit': 'named workload call',
                     'blocks': []} for name in names}
    for round_index in range(args.rounds):
        for name in names if round_index % 2 == 0 else reversed(names):
            result = core.measure(calls[name], args.trials, args.target_seconds)
            result['round'] = round_index
            results[name]['blocks'].append(result)
        print('round %d/%d complete' % (round_index + 1, args.rounds), flush=True)
    for name in names:
        if name.startswith('bezier_subdivide') or name == 'bezier_path_8_segments':
            results[name]['output_points'] = len(calls[name]())
    fingerprint = hashlib.sha256()
    for file in sorted((core.ROOT / 'gem').rglob('*.py')):
        fingerprint.update(file.relative_to(core.ROOT).as_posix().encode())
        fingerprint.update(file.read_bytes())
    try:
        commit = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=core.ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        commit = None
    data = {'schema': SCHEMA, 'recorded_at_utc': datetime.now(timezone.utc).isoformat(), 'environment': environment(),
            'source': {'commit': commit, 'core_sha256': fingerprint.hexdigest(),
                       'runner_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            'methodology': dict(METHOD, rounds=args.rounds, trials=args.trials,
                                target_seconds=args.target_seconds, setup_seconds=setup,
                                order='reverse workload order on alternating rounds',
                                seed='fixed source constants; no RNG'),
            'suite': args.suite, 'results': results,
            'unsupported': {'ray_intersections': 'supported Ray has no intersection method'}}
    if args.profiles:
        data['profiles'] = {name: core.profile(calls[name], max(1, min(1000, int(
            .05 / (results[name]['blocks'][0]['median_us'] / 1e6))))) for name in names}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
    print('%d workloads saved to %s' % (len(names), args.output))


if __name__ == '__main__':
    main()
