"""Matched base/repair measurements; setup and reference calculations are untimed."""
import argparse
import gc
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import plane, vector

BASE = '495dc3d5b4d5f8bd4969b690b3a207a9f3339023'


def historical(name):
    source = subprocess.check_output(['git', 'show', BASE + ':gem/' + name + '.py'], cwd=ROOT)
    module = types.ModuleType('baseline_' + name)
    exec(compile(source, 'baseline_' + name + '.py', 'exec'), module.__dict__)
    return module, hashlib.sha256(source).hexdigest()


def V(values):
    return vector.Vector(len(values), list(values))


def summary(samples):
    median = statistics.median(samples)
    return dict(median=median, mad=statistics.median(abs(x-median) for x in samples),
                minimum=min(samples), maximum=max(samples), trials=samples)


def ordinary_comparison(old_vector, old_plane):
    rng = random.Random(4404)
    records = {}
    for equal in (False, True):
        identical = 0
        maximum = 0.
        for _ in range(1000):
            angle = rng.uniform(0., math.pi/3.)
            incident = V([math.sin(angle), -math.cos(angle), 0.])
            normal = V([0., 1., 0.])
            eta = 1. if equal else rng.uniform(.5, 1.1)
            before = old_vector.refract(eta, incident, normal).vector
            after = vector.refract(eta, incident, normal).vector
            identical += before == after
            maximum = max(maximum, *(abs(x-y) for x,y in zip(before, after)))
        records['equal_indices' if equal else 'other_indices'] = dict(
            cases=1000, bit_identical=identical, max_absolute_difference=maximum)
    identical = 0
    maximum = 0.
    for _ in range(1000):
        points = [V([rng.uniform(-10., 10.) for axis in range(3)]) for vertex in range(4)]
        before = old_plane.Plane().bestFitNormal(points).vector
        after = plane.Plane().bestFitNormal(points).vector
        identical += before == after
        maximum = max(maximum, *(abs(x-y) for x,y in zip(before, after)))
    records['ordinary_newell'] = dict(cases=1000, bit_identical=identical,
                                     max_absolute_difference=maximum)
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--trials', type=int, default=9)
    parser.add_argument('--iterations', type=int, default=20000)
    args = parser.parse_args()
    if args.trials < 3 or args.iterations < 1:
        parser.error('use at least three trials and positive iterations')
    old_vector, vector_hash = historical('vector')
    old_plane, plane_hash = historical('plane')
    normal = V([0., 1., 0.])
    workloads = []
    for name, eta, incident in [
        ('equal_ordinary', 1., [.6, -.8, 0.]),
        ('equal_grazing', 1., [1., -1e-100, 0.]),
        ('snell_air_glass', 2./3., [.6, -.8, 0.]),
        ('snell_glass_air', 1.5, [.6, -.8, 0.]),
        ('total_internal_reflection', 1.5, [.9, -math.sqrt(.19), 0.]),
        ('critical_boundary', 1.25, [.8, -.6, 0.]),
    ]:
        incident = V(incident)
        workloads.append((name, args.iterations,
                          lambda m=old_vector, e=eta, i=incident: m.refract(e, i, normal),
                          lambda e=eta, i=incident: vector.refract(e, i, normal)))
    for count in (3, 4, 8, 64, 256):
        points = [V([3.*math.cos(2.*math.pi*i/count), 2.*math.sin(2.*math.pi*i/count),
                     math.cos(2.*math.pi*i/count) + .5*math.sin(2.*math.pi*i/count)])
                  for i in range(count)]
        old_helper, helper = old_plane.Plane(), plane.Plane()
        workloads.append(('polygon_' + str(count), max(1, args.iterations//max(1, count//4)),
                          lambda p=points, h=old_helper: h.bestFitNormal(p),
                          lambda p=points, h=helper: h.bestFitNormal(p)))
    # Neither historical grazing outputs nor polygon outputs are correctness oracles.
    result = dict(schema_version=1, base=BASE, methodology=dict(
        clock='process_time_ns', warmup_calls=100, alternating_order=True,
        gc_disabled_during_timing=True, dataset_seed=4404, trials=args.trials,
        setup_excluded=True, call_loop_and_output_allocations_included=True),
        environment=dict(python=sys.version, implementation=platform.python_implementation(),
                         platform=platform.platform(), machine=platform.machine(),
                         processor=platform.processor(), six_version=__import__('six').__version__),
        source_sha256=dict(base_vector=vector_hash, base_plane=plane_hash,
                           repaired_vector=hashlib.sha256((ROOT/'gem/vector.py').read_bytes()).hexdigest(),
                           repaired_plane=hashlib.sha256((ROOT/'gem/plane.py').read_bytes()).hexdigest()),
        ordinary_comparison=ordinary_comparison(old_vector, old_plane), workloads={})
    cpuinfo = Path('/proc/cpuinfo')
    if cpuinfo.exists():
        result['environment']['cpu_model'] = next((line.split(':', 1)[1].strip()
            for line in cpuinfo.read_text().splitlines() if line.startswith('model name')), 'unknown')
    for name, count, before, after in workloads:
        for _ in range(100):
            before()
            after()
        timings = [[], []]
        was_enabled = gc.isenabled()
        gc.disable()
        try:
            for trial in range(args.trials):
                for index in ((0, 1) if trial % 2 == 0 else (1, 0)):
                    function = (before, after)[index]
                    start = time.process_time_ns()
                    for _ in range(count):
                        function()
                    timings[index].append((time.process_time_ns()-start)/count)
        finally:
            if was_enabled:
                gc.enable()
        result['workloads'][name] = dict(calls_per_trial=count, unit='ns/call',
            before=summary(timings[0]), after=summary(timings[1]),
            paired_latency_ratio=summary([a/b for a,b in zip(timings[1], timings[0])]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print('Wrote', args.output)


if __name__ == '__main__':
    main()
