"""Descriptive paired-block bootstrap for saved Vector/Quaternion runs."""
import argparse
import json
from pathlib import Path
import random
import statistics


def summarize(blocks):
    ratios=[b['before']['median_us']/b['after']['median_us'] for b in blocks]
    rng=random.Random(303)
    bootstrap=sorted(statistics.median(rng.choices(ratios,k=len(ratios))) for _ in range(10000))
    return {'before_median_us':statistics.median(b['before']['median_us'] for b in blocks),
            'after_median_us':statistics.median(b['after']['median_us'] for b in blocks),
            'paired_median_speedup':statistics.median(ratios),'paired_speedup_range':[min(ratios),max(ratios)],
            'descriptive_bootstrap_95_interval':[bootstrap[249],bootstrap[9749]],
            'before_median_relative_mad_percent':statistics.median(b['before']['mad_us']/b['before']['median_us']*100 for b in blocks),
            'after_median_relative_mad_percent':statistics.median(b['after']['mad_us']/b['after']['median_us']*100 for b in blocks),
            'round_speedups':ratios,'all_blocks_faster':all(r>1 for r in ratios)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('runs',nargs='+',type=Path)
    parser.add_argument('--follow-up',type=Path)
    args=parser.parse_args()
    runs=[json.loads(path.read_text()) for path in args.runs]
    assert len(runs)>=2
    assert all(r['baseline_commit']==runs[0]['baseline_commit'] and r['current_sha256']==runs[0]['current_sha256'] for r in runs)
    cases={name:summarize([block[name] for run in runs for block in run['rounds']]) for name in runs[0]['summary']}
    result={'runs':len(runs),'blocks_per_case':sum(len(r['rounds']) for r in runs),
            'method':'10000 seeded paired-round resamples; percentile interval for median speedup',
            'limitations':'descriptive conditional interval; block independence and stationary shared-host load not guaranteed; not a population guarantee',
            'source_sha256':runs[0]['current_sha256'],'cases':cases}
    if args.follow_up:
        follow=json.loads(args.follow_up.read_text())
        assert follow['baseline_commit']==runs[0]['baseline_commit'] and follow['current_sha256']==runs[0]['current_sha256']
        result['follow_up']={'source':args.follow_up.name,'blocks':len(follow['rounds']),
                            'target_seconds':follow['methodology']['target_seconds'],
                            'reason':'investigate unchanged-control shifts and inconclusive squad4 result',
                            'cases':{name:summarize([b[name] for b in follow['rounds']]) for name in follow['summary']}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(f"{len(cases)} cases summarized from {result['blocks_per_case']} primary blocks")


if __name__=='__main__':main()
