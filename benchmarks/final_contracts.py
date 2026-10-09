"""Paired Phase 2G equality/clamp measurements; no unrelated optimization."""
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import vector
from benchmarks.core_baseline import measure

BASE = 'c1e3fe745668b205dbd7e2247350fbb6b4f0962c'


def main():
    namespace = {'__name__': 'baseline_vector'}
    text = subprocess.check_output(['git', 'show', BASE + ':gem/vector.py'], cwd=ROOT, text=True)
    exec(compile(text, 'baseline_vector.py', 'exec'), namespace)
    cases = {}
    for size in (2, 3, 4):
        values = [float(i + 1) for i in range(size)]
        changed = values[:]; changed[-1] += .25
        for label, operand in [('equal', values), ('different_last', changed)]:
            a, b = namespace['Vector'](size, values[:]), namespace['Vector'](size, operand[:])
            c, d = vector.Vector(size, values[:]), vector.Vector(size, operand[:])
            cases[f'equality_{size}_{label}'] = (lambda a=a,b=b:a == b, lambda c=c,d=d:c == d)
            cases[f'inequality_{size}_{label}'] = (lambda a=a,b=b:a != b, lambda c=c,d=d:c != d)
        inputs = [-2, 2, 10, -8][:size]; lower, upper = [0]*size, [5]*size
        old, new = namespace['clamp'], vector.clamp
        # Both variants receive an identical fresh list on every call: resetting
        # outside timeit would let the historical mutation alter subsequent inputs.
        cases[f'clamp_{size}'] = (lambda old=old,size=size,v=inputs,lo=lower,hi=upper:old(size,v[:],lo,hi),
                                lambda new=new,size=size,v=inputs,lo=lower,hi=upper:new(size,v[:],lo,hi))
    rounds = []
    for round_index in range(3):
        results = {}
        for name, pair in cases.items():
            order = (0, 1) if round_index % 2 == 0 else (1, 0)
            times = {}
            for index in order:times['before' if index == 0 else 'after'] = measure(pair[index], 7, .05)
            times['after_before_ratio'] = times['after']['median_us']/times['before']['median_us']
            results[name] = times
        rounds.append(results)
    summary = {name: {'before_median_us': statistics.median(r[name]['before']['median_us'] for r in rounds),
                      'after_median_us': statistics.median(r[name]['after']['median_us'] for r in rounds),
                      'paired_median_ratio': statistics.median(r[name]['after_before_ratio'] for r in rounds),
                      'paired_ratio_range': [min(r[name]['after_before_ratio'] for r in rounds),
                                             max(r[name]['after_before_ratio'] for r in rounds)]}
               for name in cases}
    output = {'baseline_commit': BASE, 'python': sys.version, 'platform': platform.platform(),
              'methodology': '3 rounds; alternating before/after order; 7 trials; .05s calibration; same process',
              'input_ranges': 'Vector2/3/4 values 1..4; last component +.25; clamp inputs -8..10, bounds 0..5',
              'clamp_input_reset': 'fresh value-list copy included equally in both variants; API result allocation included',
              'new_domains': 'empty/mixed-size comparisons are correctness changes, not comparable old boolean operations',
              'summary': summary, 'rounds': rounds}
    dest = Path(sys.argv[1]); dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(output, indent=2)+'\n')
    print(f'{len(cases)} paired cases saved to {dest}')


if __name__ == '__main__':
    main()
