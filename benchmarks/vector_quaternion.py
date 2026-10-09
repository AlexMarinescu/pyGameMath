"""Paired Vector/Quaternion comparison against merged Phase 3B master."""
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

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import vector, quaternion
from benchmarks import core_baseline as core

BASE = '86dcd405da75c10da976ac67c0b02d774f066920'


def baseline():
    modules, hashes = [], {}
    for name in ('vector', 'quaternion'):
        source = subprocess.check_output(['git', 'show', BASE+':gem/'+name+'.py'], cwd=ROOT, text=True)
        mod = types.ModuleType('baseline_'+name)
        exec(compile(source,'baseline_'+name+'.py','exec'),vars(mod))
        modules.append(mod)
        hashes[name] = hashlib.sha256(source.encode()).hexdigest()
    modules[1].vector = modules[0]
    return modules, hashes


def cases(v, q):
    # Compile the unchanged Phase 3A section with independent module globals.
    source = (ROOT/'benchmarks/core_baseline.py').read_text()
    source = source[source.index('def prepare():'):source.index('    translation=matrix.Matrix(4)')]
    namespace = dict(vars(core),vector=v,quaternion=q,
                     V=lambda *values:v.Vector(len(values),list(values)),
                     V3=v.Vector(3,[1.25,-2.5,3.75]), V3B=v.Vector(3,[-.5,2,1]))
    exec(compile(source+'    return calls, details\n','phase3a_vector_quaternion.py','exec'),namespace)
    all_calls, details = namespace['prepare']()
    calls = {k: fn for k, fn in all_calls.items() if k.startswith(('vector','quaternion'))}
    for n in (2,3,4):
        values = [(-1.)**i*(i+.25) for i in range(n)]
        calls[f'vector{n}_raw_normalize'] = lambda n=n,values=values:v.normalize(n,values)
    a = q.Quaternion([.5,.5,.5,.5]); b = q.Quaternion([.5,-.5,.5,-.5])
    calls['quaternion_raw_multiply'] = lambda:q.quat_mul_quat(a.data,b.data)
    calls['quaternion_raw_conjugate'] = lambda:q.quat_conjugate(a.data)
    calls['quaternion_conjugate'] = a.conjugate
    calls['quaternion_raw_slerp'] = lambda:q.quat_slerp(a,b,.37)
    values=[-2.,.5,10.];lo=[0.]*3;hi=[5.]*3
    calls['vector3_clamp'] = lambda:v.clamp(3,values,lo,hi)
    calls['vector3_in_place_normalize'] = lambda:v.Vector(3,list(values)).i_normalize()
    rows=[[1.,0,0,0],[0,1.,0,0],[0,0,1.,0],[2.,-3.,.5,1.]]
    calls['vector3_affine_transform'] = lambda:v.transform(3,values,rows)
    calls['vector3_in_place_transform'] = lambda:v.Vector(3,list(values)).i_transform(values,rows)
    calls['vector3_normalize_zero'] = lambda:v.normalize(3,[0.,0.,0.])
    calls['quaternion_normalize_zero'] = lambda:q.quat_normalize([0.,0.,0.,0.])
    return calls, details


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rounds',type=int,default=3)
    parser.add_argument('--trials',type=int,default=7)
    parser.add_argument('--target-seconds',type=float,default=.02)
    parser.add_argument('--case',action='append',dest='selected',help='repeat to restrict a follow-up to named cases')
    args=parser.parse_args()
    if args.rounds<3 or args.trials<5 or not math.isfinite(args.target_seconds) or args.target_seconds<=0:
        parser.error('require >=3 rounds, >=5 trials and positive finite target')
    modules, hashes=baseline()
    before,details=cases(*modules);after,_=cases(vector,quaternion)
    if args.selected:
        unknown=set(args.selected)-set(before)
        if unknown:parser.error('unknown cases: '+', '.join(sorted(unknown)))
        before={name:before[name] for name in args.selected}
        after={name:after[name] for name in args.selected}
    rounds=[]
    for index in range(args.rounds):
        rows={}
        for name in before:
            pair={}
            for label in (('before','after') if index%2==0 else ('after','before')):
                pair[label]=core.measure((before if label=='before' else after)[name],args.trials,args.target_seconds)
            rows[name]=pair
        rounds.append(rows)
        print(f'Completed round {index+1}/{args.rounds}',flush=True)
    summary={}
    for name in before:
        ratios=[r[name]['before']['median_us']/r[name]['after']['median_us'] for r in rounds]
        summary[name]={'before_median_us':statistics.median(r[name]['before']['median_us'] for r in rounds),
                       'after_median_us':statistics.median(r[name]['after']['median_us'] for r in rounds),
                       'paired_median_speedup':statistics.median(ratios),'paired_speedup_range':[min(ratios),max(ratios)]}
    selected=('vector3_normalize','vector3_cross','quaternion_rotate_vector','quaternion_slerp',
              'quaternion_squad','quaternion_squad4','quaternion_multiply','quaternion_to_matrix')
    profiles={name:{label:core.profile(fn[name],1000) for label,fn in [('before',before),('after',after)]} for name in selected if name in before}
    import six
    cpu=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),'unavailable')
    data={'baseline_commit':BASE,'baseline_sha256':hashes,
          'current_sha256':{n:hashlib.sha256((ROOT/'gem'/f'{n}.py').read_bytes()).hexdigest() for n in ('vector','quaternion')},
          'environment':{'python':sys.version,'implementation':platform.python_implementation(),'platform':platform.platform(),'cpu':cpu,'six':six.__version__},
          'methodology':{'rounds':args.rounds,'trials':args.trials,'target_seconds':args.target_seconds,
                         'order':'alternate before/after per round','timer':'perf_counter, timeit disables GC',
                         'setup':'Phase 3A deterministic input generation excluded; result allocation included',
                         'in_place':'fresh receiver construction/copy included to avoid input drift',
                         'profiles':'1000 calls, separately instrumented; traced peak per call is not RSS'},
          'datasets':details,'rounds':rounds,'summary':summary,'profiles':profiles}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2)+'\n')
    print(f'{len(before)} paired cases saved to {args.output}')


if __name__=='__main__':main()
