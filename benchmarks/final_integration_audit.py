"""Measure the specific whole-path validation hypothesis; no timing gate.

Run: python benchmarks/final_integration_audit.py --output /tmp/path-cost.json
Setup is excluded. Straight cubic segments require no adaptive subdivision.
"""
import argparse
import cProfile
import gc
import json
import platform
from pathlib import Path
import random
import statistics
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem.bezier import BezierPath
from gem.vector import Vector


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    workloads = {}
    for segments, count in ((8, 16), (32, 8), (128, 2), (512, 1)):
        path = BezierPath()
        path.setControlPoints([Vector(2, [float(i), 0.]) for i in range(3*segments+1)])
        workloads[segments] = (path, count)
        actual = [point for curve in path.getDrawingPoints() for point in curve]
        assert [point.vector for point in actual] == [[float(3*i), 0.] for i in range(segments+1)]
        assert all(point is not control and point.vector is not control.vector
                   for point in actual for control in path.controlPoints)
    measurements = {n: [] for n in workloads}
    rng = random.Random(475005)
    for trial in range(7):
        order = list(workloads)
        rng.shuffle(order)
        for segments in order:
            path, count = workloads[segments]
            was_enabled = gc.isenabled()
            gc.disable()
            try:
                start = time.process_time_ns()
                for _ in range(count):
                    path.getDrawingPoints()
                elapsed = time.process_time_ns()-start
            finally:
                if was_enabled:
                    gc.enable()
            measurements[segments].append(elapsed/count)
    rows = []
    for segments, (path, count) in workloads.items():
        profiler = cProfile.Profile()
        profiler.runcall(path.getDrawingPoints)
        entries = [entry for entry in profiler.getstats()
                   if hasattr(entry.code, 'co_name') and entry.code.co_name == '_coordinates']
        assert len(entries) == 1 and entries[0].callcount == segments+1
        values = measurements[segments]
        median = statistics.median(values)
        rows.append({'segments': segments, 'controls': 3*segments+1, 'output_points': segments+1,
                     'calls_per_trial': count, 'trials_ns_per_call': values, 'median_ns': median,
                     'mad_ns': statistics.median(abs(x-median) for x in values),
                     'min_ns': min(values), 'max_ns': max(values),
                     'profile_coordinates_calls': entries[0].callcount,
                     'profile_coordinates_seconds': entries[0].totaltime,
                     'profile_total_seconds': max(e.totaltime for e in profiler.getstats()),
                     'validated_control_visits': (segments+1)*(3*segments+1)})
    result = {'schema': 1, 'hypothesis': 'whole-path validation repeated for each straight segment',
              'environment': {'python': sys.version, 'platform': platform.platform(),
                              'machine': platform.machine()},
              'method': {'clock': 'process_time_ns', 'trials': 7, 'seed': 475005,
                         'gc_disabled_during_timing': True, 'warmup_calls': 1,
                         'setup_timed': False, 'profile_separate_from_timing': True,
                         'control_range': 'X=0..3*segments; Y=0; Vector2',
                         'subdivision': 'straight segments, endpoints only'}, 'workloads': rows}
    args.output.write_text(json.dumps(result, indent=2)+'\n')
    for row in rows:
        print('{} segments: {:.3f} ms, MAD {:.3f} ms; {} control visits'.format(
            row['segments'], row['median_ns']/1e6, row['mad_ns']/1e6, row['validated_control_visits']))


if __name__ == '__main__':
    main()
