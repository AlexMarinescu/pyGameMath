"""Paired power/log timings and unguarded-squaring drift for Phase 4G-1R.

Run from a checkout; git supplies the immutable pre-repair implementation.
Setup/module loading is excluded. No installed-runtime dependencies are added.
"""
import argparse
from decimal import Decimal, localcontext
import gc
import hashlib
import json
import math
import platform
from pathlib import Path
import random
import statistics
import struct
import subprocess
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import quaternion  # noqa: E402

BASE = 'f540db6cc5d48ec9fdd60d0972f26f636c66ab19'


def previous_module(base):
    source = subprocess.check_output(['git', 'show', base+':gem/quaternion.py'],
                                     cwd=str(ROOT), text=True)
    module = types.ModuleType('quaternion_before_repairs')
    exec(compile(source, base+':gem/quaternion.py', 'exec'), module.__dict__)
    return module


def verify_ordinary_compatibility(previous):
    """Check actual baseline code, not just the equivalent reference formula."""
    comparisons = 0
    for seed in range(16):
        rng = random.Random(471000+seed)
        axis = [rng.uniform(-2, 2) for _ in range(3)]
        length = math.sqrt(sum(x*x for x in axis))
        angle = rng.uniform(0.1, 2.9)
        data = [math.cos(angle)] + [x/length*math.sin(angle) for x in axis]
        for scale in [1.0, 2.0]:
            values = [x*scale for x in data]
            before, after = previous.Quaternion(values[:]), quaternion.Quaternion(values[:])
            for exponent in [-16, -1, -0.5, 0, 0.5, 1, 2, 16]:
                a, b = before.pow(exponent).data, after.pow(exponent).data
                assert [struct.pack('!d', x) for x in a] == [struct.pack('!d', x) for x in b]
                comparisons += 1
            assert [struct.pack('!d', x) for x in before.log()] == [struct.pack('!d', x) for x in after.log()]
            comparisons += 1
    return comparisons


def unguarded_power(data, exponent):
    """Investigation only; intentionally not the supported core dispatch."""
    result, base = [1.0, 0.0, 0.0, 0.0], list(data)
    while exponent:
        if exponent & 1:
            result = quaternion.quat_mul_quat(result, base)
        exponent >>= 1
        if exponent:
            base = quaternion.quat_mul_quat(base, base)
    return result


def drift_measurements():
    rows = []
    for name, data in [
        ('basis_exact', [0.0, 1.0, 0.0, 0.0]),
        ('halves_exact', [0.5]*4),
        ('general_angle_0.3', [math.cos(0.3), math.sin(0.3), 0.0, 0.0]),
        ('rounded_eighth_turn', [-math.sqrt(0.5), math.sqrt(0.5), 0.0, 0.0]),
    ]:
        with localcontext() as context:
            context.prec = 100
            exact_squared = sum(Decimal.from_float(x)**2 for x in data)
            offset = str(exact_squared.sqrt()-1)
        for exponent in [32, 1024, 10**6, 10**16]:
            candidate = unguarded_power(data, exponent)
            current = quaternion.Quaternion(list(data)).pow(exponent).data
            rows.append({'name': name, 'input': data, 'exponent': exponent,
                         'represented_input_norm_minus_one_decimal': offset,
                         'unguarded_output': candidate,
                         'unguarded_norm': math.hypot(*candidate),
                         'selected_output': current,
                         'selected_norm': math.hypot(*current)})
    return rows


def median_absolute_deviation(values):
    median = statistics.median(values)
    return statistics.median(abs(x-median) for x in values)


def summarize(values):
    return {'median_ns': statistics.median(values),
            'mad_ns': median_absolute_deviation(values),
            'min_ns': min(values), 'max_ns': max(values)}


def latency(call, count):
    start = time.process_time_ns()
    for _ in range(count):
        call()
    return (time.process_time_ns()-start)/count


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base', default=BASE)
    parser.add_argument('--trials', type=int, default=9)
    parser.add_argument('--iterations', type=int, default=30000)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.trials < 3 or args.iterations < 100:
        parser.error('Use at least three trials and 100 iterations')
    previous = previous_module(args.base)
    compatibility_comparisons = verify_ordinary_compatibility(previous)
    general = [math.cos(0.4)] + [x*math.sin(0.4) for x in [1/3, -2/3, 2/3]]
    tiny = math.ldexp(1.0, -1074)
    definitions = [
        ('pow_general_8', 'pow', general, 8),
        ('pow_general_half', 'pow', general, 0.5),
        ('pow_copy_1', 'pow', general, 1),
        ('pow_nonunit_2', 'pow', [2.0, 1.0, 0.0, 0.0], 2),
        ('pow_cyclic_8', 'pow', [0.5]*4, 8),
        ('pow_cyclic_large', 'pow', [0.5]*4, 10**16),
        ('pow_subnormal_half', 'pow', [-1.0, tiny, tiny, 0.0], 0.5),
        ('log_general', 'log', general, None),
        ('log_subnormal', 'log', [-1.0, tiny, tiny, 0.0], None),
    ]
    workloads = []
    for name, method, data, exponent in definitions:
        before, after = previous.Quaternion(data[:]), quaternion.Quaternion(data[:])
        before_call, after_call = getattr(before, method), getattr(after, method)
        if exponent is not None:
            before_call = lambda call=before_call, e=exponent: call(e)
            after_call = lambda call=after_call, e=exponent: call(e)
        # Raw function versus wrapper overhead is measured separately for the
        # ordinary paths; the two runs use the same input/data-generation costs.
        entry = {'name': name, 'method': method, 'input': data,
                 'exponent': exponent, 'before': [], 'after': [], 'pairs': [],
                 'operations_per_trial': max(100, args.iterations//15) if name == 'pow_cyclic_large' else args.iterations,
                 'before_call': before_call, 'after_call': after_call}
        workloads.append(entry)
        if name in ['pow_general_half', 'log_general']:
            old_function = previous.quat_pow if method == 'pow' else previous.quat_log
            new_function = quaternion.quat_pow if method == 'pow' else quaternion.quat_log
            raw = {key: value for key, value in entry.items() if key not in ['before_call', 'after_call']}
            raw.update(name=name+'_raw', before=[], after=[], pairs=[])
            if method == 'pow':
                raw['before_call'] = lambda f=old_function, q=before, e=exponent: f(q, e)
                raw['after_call'] = lambda f=new_function, q=after, e=exponent: f(q, e)
            else:
                raw['before_call'] = lambda f=old_function, q=before: f(q)
                raw['after_call'] = lambda f=new_function, q=after: f(q)
            workloads.append(raw)
    was_enabled = gc.isenabled()
    gc.disable()
    try:
        for work in workloads:
            for _ in range(200):
                work['before_call']()
                work['after_call']()
        for trial in range(args.trials):
            order = list(workloads)
            random.Random(471001+trial).shuffle(order)
            for work in order:
                directions = ['before', 'after'] if trial % 2 == 0 else ['after', 'before']
                pair = {}
                for side in directions:
                    elapsed = latency(work[side+'_call'], work['operations_per_trial'])
                    work[side].append(elapsed)
                    pair[side] = elapsed
                work['pairs'].append(pair['after']/pair['before'])
    finally:
        if was_enabled:
            gc.enable()
    for work in workloads:
        work['before_summary'] = summarize(work['before'])
        work['after_summary'] = summarize(work['after'])
        work['paired_after_before'] = {'median': statistics.median(work['pairs']),
                                     'mad': median_absolute_deviation(work['pairs']),
                                     'min': min(work['pairs']), 'max': max(work['pairs'])}
        del work['before_call'], work['after_call']
    output = {'schema': 1, 'baseline_commit': args.base,
              'ordinary_compatibility_comparisons': compatibility_comparisons,
              'measured_quaternion_sha256': hashlib.sha256((ROOT/'gem/quaternion.py').read_bytes()).hexdigest(),
              'environment': {'python': sys.version, 'implementation': platform.python_implementation(),
                              'platform': platform.platform(), 'machine': platform.machine()},
              'methodology': {'clock': 'process_time_ns', 'trials': args.trials,
                              'setup_excluded': True, 'warmup_calls_per_side': 200,
                              'gc_disabled_during_timings': True, 'loop_overhead_subtracted': False,
                              'order': 'alternate before/after; workload shuffle seed 471001+trial'},
              'workloads': workloads, 'unguarded_squaring_drift': drift_measurements()}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False)+'\n')
    for work in workloads:
        print('{}: {:.3f} -> {:.3f} us, paired ratio {:.3f}'.format(
            work['name'], work['before_summary']['median_ns']/1000,
            work['after_summary']['median_ns']/1000,
            work['paired_after_before']['median']))


if __name__ == '__main__':
    main()
