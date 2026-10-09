"""Interleaved post-2G inverse comparison; retains numerical stabilization."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import matrix
from benchmarks.core_baseline import measure, profile

BASE = '89cd4b97784d32625de65b8f465cfcbf6e102943'


def datasets(size):
    rng = random.Random(307 + size)
    ordinary = [[[3.*(i==j)+.13*(i+1)/(j+1) for j in range(size)] for i in range(size)]]
    for _ in range(3):
        ordinary.append([[(4 if i==j else 0) + rng.uniform(-.5, .5) for j in range(size)] for i in range(size)])
    yield 'ordinary', ordinary
    for scale in (1e-300, 1e300):
        yield f'scale_{scale:g}', [[[v*scale for v in row] for row in rows] for rows in ordinary]
    exponents = [-300, 300, -200, 200][:size]
    import math
    # Bidiagonal mixed-row dataset; every required identity product representable.
    yield 'mixed_rows', [[[math.ldexp(float(2 if i==j else 1 if j==i+1 else 0), exponents[i])
                           for j in range(size)] for i in range(size)]]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--rounds', type=int, default=3)
    parser.add_argument('--trials', type=int, default=7)
    parser.add_argument('--target-seconds', type=float, default=.04)
    args = parser.parse_args()
    if args.rounds<3 or args.trials<5 or not 0 < args.target_seconds < float('inf'):parser.error('require >=3 rounds, >=5 trials and finite positive target')
    text = subprocess.check_output(['git', 'show', BASE+':gem/matrix.py'], cwd=ROOT, text=True)
    ns = {'__name__': 'baseline_matrix'};exec(compile(text,'baseline_matrix.py','exec'),ns)
    cases, metadata = {}, {}
    for size in (3, 4):
        for regime, rows_list in datasets(size):
            for mode in ('raw', 'returning', 'in_place'):
                pair = []
                for namespace, cls in [(ns, ns['Matrix']), (vars(matrix), matrix.Matrix)]:
                    if mode == 'raw':
                        fn = namespace['inverse'+str(size)]
                        call = lambda fn=fn,rows_list=rows_list:[fn(rows) for rows in rows_list]
                    elif mode == 'returning':
                        objects = [cls(size, rows) for rows in rows_list]
                        call = lambda objects=objects:[obj.inverse() for obj in objects]
                    else:
                        # Include receiver construction equally to avoid inverse-of-inverse drift.
                        call = lambda cls=cls,size=size,rows_list=rows_list:[cls(size,rows).i_inverse() for rows in rows_list]
                    pair.append(call)
                name = f'matrix{size}_{regime}_{mode}'
                cases[name] = pair
                metadata[name] = {'size': size, 'regime': regime, 'mode': mode, 'inverses_per_call':len(rows_list)}
    rounds = []
    for iteration in range(args.rounds):
        result = {}
        for name, pair in cases.items():
            timing = {}
            for index in ((0,1) if iteration%2==0 else (1,0)):
                value = measure(pair[index],args.trials,args.target_seconds)
                count = metadata[name]['inverses_per_call']
                for field in ('median_us','min_us','max_us','mad_us'):value[field]/=count
                value['seconds_per_operation'] = [v/count for v in value['seconds_per_operation']]
                value['inversions_per_trial'] = value['operations_per_trial']*count
                timing['before' if index==0 else 'after'] = value
            result[name] = timing
        rounds.append(result)
    summaries = {}
    for name in cases:
        ratios = [r[name]['before']['median_us']/r[name]['after']['median_us'] for r in rounds]
        summaries[name] = {**metadata[name],
            'before_median_us':statistics.median(r[name]['before']['median_us'] for r in rounds),
            'after_median_us':statistics.median(r[name]['after']['median_us'] for r in rounds),
            'paired_median_speedup':statistics.median(ratios),'paired_speedup_range':[min(ratios),max(ratios)],
            'before_round_medians_us':[r[name]['before']['median_us'] for r in rounds],
            'after_round_medians_us':[r[name]['after']['median_us'] for r in rounds]}
    profiles = {}
    for name in ('matrix3_ordinary_raw','matrix4_ordinary_raw','matrix4_ordinary_returning',
                 'matrix4_scale_1e+300_raw','matrix3_mixed_rows_raw'):
        profiles[name] = {label:profile(cases[name][i],100) for i,label in enumerate(('before','after'))}
    import six
    cpu = 'unavailable'
    if Path('/proc/cpuinfo').exists():
        cpu = next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),cpu)
    result = {'baseline_commit':BASE,'baseline_source_sha256':hashlib.sha256(text.encode()).hexdigest(),
              'current_source_sha256':hashlib.sha256((ROOT/'gem/matrix.py').read_bytes()).hexdigest(),
              'environment':{'python':sys.version,'implementation':platform.python_implementation(),
                             'platform':platform.platform(),'cpu':cpu,'six':six.__version__},
              'methodology':{'rounds':args.rounds,'trials':args.trials,'target_seconds':args.target_seconds,
                             'order':'alternate before/after each round','timer':'perf_counter, timeit GC disabled',
                             'timings':'per inverse, normalized by batch size; driver/list allocation included',
                             'in_place':'fresh receiver construction included equally',
                             'setup':'data generation and returning-receiver construction excluded',
                             'profiles':'100 batches; instrumented times not latency; peak is per batch'},
              'summary':summaries,'rounds':rounds,'profiles':profiles}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(f'{len(cases)} paired inverse cases saved to {args.output}')


if __name__ == '__main__':main()
