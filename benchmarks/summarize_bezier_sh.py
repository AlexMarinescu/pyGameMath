"""Summarize paired Bezier/SH timing blocks with descriptive intervals."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from benchmarks.summarize_vector_quaternion import summarize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path)
    parser.add_argument('runs',nargs='+',type=Path)
    args=parser.parse_args();runs=[json.loads(path.read_text()) for path in args.runs]
    assert len(runs)>=2
    assert all(r['baseline_commit']==runs[0]['baseline_commit'] and r['current_sha256']==runs[0]['current_sha256']
               and r['sample_counts']==runs[0]['sample_counts'] for r in runs)
    result={'runs':len(runs),'blocks_per_case':sum(len(r['rounds']) for r in runs),
            'baseline_commit':runs[0]['baseline_commit'],'source_sha256':runs[0]['current_sha256'],
            'method':'10000 seeded paired-round resamples; percentile interval for median speedup',
            'limitations':'descriptive conditional interval; shared-host block independence and stationarity not assured; no population guarantee',
            'sample_counts':runs[0]['sample_counts'],
            'cases':{name:summarize([block[name] for run in runs for block in run['rounds']]) for name in runs[0]['summary']}}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(f"{len(result['cases'])} cases summarized from {result['blocks_per_case']} paired blocks")


if __name__=='__main__':main()
