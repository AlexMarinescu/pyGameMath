"""Descriptive paired-block bootstrap for saved inverse benchmark runs."""
import json
from pathlib import Path
import random
import statistics
import sys


def main():
    output = Path(sys.argv[1])
    runs = [json.loads(Path(path).read_text()) for path in sys.argv[2:]]
    assert len(runs) >= 2
    assert all(run['baseline_commit']==runs[0]['baseline_commit'] and
               run['current_source_sha256']==runs[0]['current_source_sha256'] for run in runs)
    cases = {}
    for name in runs[0]['summary']:
        blocks = [block[name] for run in runs for block in run['rounds']]
        ratios = [b['before']['median_us']/b['after']['median_us'] for b in blocks]
        rng = random.Random(303)
        bootstrap = sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(10000))
        cases[name] = {'before_median_us':statistics.median(b['before']['median_us'] for b in blocks),
                       'after_median_us':statistics.median(b['after']['median_us'] for b in blocks),
                       'paired_median_speedup':statistics.median(ratios),'paired_speedup_range':[min(ratios),max(ratios)],
                       'descriptive_bootstrap_95_interval':[bootstrap[249],bootstrap[9749]],
                       'before_median_relative_mad_percent':statistics.median(b['before']['mad_us']/b['before']['median_us']*100 for b in blocks),
                       'after_median_relative_mad_percent':statistics.median(b['after']['mad_us']/b['after']['median_us']*100 for b in blocks),
                       'round_speedups':ratios,'all_blocks_faster':all(r>1 for r in ratios)}
    result = {'runs':len(runs),'blocks_per_case':sum(len(r['rounds']) for r in runs),
              'method':'10000 seeded paired-round resamples; percentile interval for median speedup',
              'limitations':'descriptive conditional interval; block independence and stationary shared-host load not guaranteed; not a population guarantee',
              'source_sha256':runs[0]['current_source_sha256'],'cases':cases}
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(f"{len(cases)} cases summarized from {result['blocks_per_case']} blocks")


if __name__=='__main__':main()
