"""Paired Bezier/SH measurements against merged Phase 3C master."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import types

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from gem import bezier, spherical_harmonics as sh
from benchmarks import core_baseline as core

BASE='09fb57a31fd1b62fad027822cb4ab1eba52ca158'


def baseline():
    modules=[];hashes={}
    for name in ('bezier','spherical_harmonics'):
        source=subprocess.check_output(['git','show',BASE+':gem/'+name+'.py'],cwd=ROOT,text=True)
        mod=types.ModuleType('baseline_'+name)
        exec(compile(source,'baseline_'+name+'.py','exec'),vars(mod))
        modules.append(mod);hashes[name]=hashlib.sha256(source.encode()).hexdigest()
    return modules,hashes


def cases(b,s):
    source=(ROOT/'benchmarks/core_baseline.py').read_text()
    source=source[source.index('def V('):source.index('\n\nV3,V3B=')]
    namespace=dict(vars(core),bezier=b,sh=s)
    exec(compile(source,'phase3a_bezier_sh.py','exec'),namespace)
    all_calls,details=namespace['prepare']()
    calls={name:fn for name,fn in all_calls.items() if name.startswith(('bezier','sh_'))}
    polygon=[(0.,0.,0.),(.2,2.,.5),(.8,-2.,-.5),(1.,0.,0.)]
    calls['bezier_raw_split']=lambda:b._split(polygon)
    calls['bezier_raw_subdivide']=lambda:b._subdivide(polygon,.01)
    calls['bezier_raw_flatness']=lambda:max(b._chord_distance(p,polygon[0],polygon[-1]) for p in polygon[1:-1]) if not hasattr(b,'_flatness') else b._flatness(polygon)
    for bands in (3,5,9):
        calls[f'sh_raw_basis_b{bands}']=lambda bands=bands:s._basis(bands,.73,1.27)
    # Angular projection measures basis generation plus compensated accumulation.
    hdr=[[[2.+(j+.5)/32, 1.+(i+.5)/16, .25] for j in range(32)] for i in range(16)]
    calls['sh_angular_probe_b3']=lambda:s.project_angular_probe(hdr)
    return calls,details


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rounds',type=int,default=3)
    parser.add_argument('--trials',type=int,default=7)
    parser.add_argument('--target-seconds',type=float,default=.02)
    parser.add_argument('--case',action='append',dest='selected')
    args=parser.parse_args()
    if args.rounds<3 or args.trials<5 or not math.isfinite(args.target_seconds) or args.target_seconds<=0:
        parser.error('require >=3 rounds, >=5 trials and positive finite target')
    modules,hashes=baseline();before,details=cases(*modules);after,_=cases(bezier,sh)
    if args.selected:
        if set(args.selected)-set(before):parser.error('unknown cases')
        before={n:before[n] for n in args.selected};after={n:after[n] for n in args.selected}
    counts={}
    for name in before:
        if name.startswith('bezier_subdivide') or name=='bezier_raw_subdivide':
            pair=[len(before[name]()),len(after[name]())]
            assert pair[0]==pair[1],name
            counts[name]={'before':pair[0],'after':pair[1]}
        elif name=='bezier_path_8_segments':
            pair=[list(map(len,before[name]())),list(map(len,after[name]()))]
            assert pair[0]==pair[1],name
            counts[name]={'before':pair[0],'after':pair[1]}
    rounds=[]
    for index in range(args.rounds):
        rows={}
        for name in before:
            pair={}
            for label in (('before','after') if index%2==0 else ('after','before')):
                pair[label]=core.measure((before if label=='before' else after)[name],args.trials,args.target_seconds)
            rows[name]=pair
        rounds.append(rows);print(f'Completed round {index+1}/{args.rounds}',flush=True)
    summary={}
    for name in before:
        ratios=[r[name]['before']['median_us']/r[name]['after']['median_us'] for r in rounds]
        summary[name]={'before_median_us':statistics.median(r[name]['before']['median_us'] for r in rounds),
                       'after_median_us':statistics.median(r[name]['after']['median_us'] for r in rounds),
                       'paired_median_speedup':statistics.median(ratios),'paired_speedup_range':[min(ratios),max(ratios)]}
    selected=('bezier_cubic','bezier_subdivide_0.0001','bezier_path_8_segments','sh_reconstruct_l2','sh_rotate_l2_rgb','sh_angular_probe_b3','sh_project_n256_b3')
    profiles={name:{label:core.profile(fn[name],5 if 'subdivide' in name or 'path' in name or 'probe' in name else 100)
                    for label,fn in [('before',before),('after',after)]} for name in selected if name in before}
    import six
    cpu=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),'unavailable')
    data={'baseline_commit':BASE,'baseline_sha256':hashes,
          'current_sha256':{n:hashlib.sha256((ROOT/'gem'/f'{n}.py').read_bytes()).hexdigest() for n in hashes},
          'environment':{'python':sys.version,'implementation':platform.python_implementation(),'platform':platform.platform(),'cpu':cpu,'six':six.__version__},
          'methodology':{'rounds':args.rounds,'trials':args.trials,'target_seconds':args.target_seconds,
                         'order':'alternate before/after per round','timer':'perf_counter, timeit disables GC',
                         'setup':'Phase 3A input generation excluded; returned result allocation included',
                         'basis_cache':'warm execution: calibration/warmup excluded; bounded 16-layout cache',
                         'profiles':'separate instrumented calls; traced peak per call is not RSS'},
          'datasets':details,'sample_counts':counts,'rounds':rounds,'summary':summary,'profiles':profiles}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n');print(f'{len(before)} paired cases saved to {args.output}')


if __name__=='__main__':main()
