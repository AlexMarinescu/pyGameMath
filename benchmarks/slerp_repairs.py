"""Matched SLERP/SQUAD4 repair measurements against immutable master source.

python benchmarks/slerp_repairs.py --output /tmp/slerp-performance.json
Setup and reference comparisons are excluded; no timing assertion or runtime dependency.
"""
import argparse
import gc
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import statistics
import struct
import subprocess
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[1]
BASE = '4e48dc68cc23d47abd7501897c6dfbb0c05156c3'
sys.path.insert(0,str(ROOT))
from gem import quaternion as q


def old_module():
    source = subprocess.check_output(['git','show',BASE+':gem/quaternion.py'],cwd=ROOT,text=True)
    module = types.ModuleType('slerp_before_repair')
    exec(compile(source,BASE+':gem/quaternion.py','exec'),module.__dict__)
    return module,source


def packed(values):
    return [struct.pack('!d',x) for x in values]


def compatibility(previous):
    interior, endpoints, changed_endpoints, legacy = 0,0,0,0
    largest_endpoint_difference = 0.
    for seed in range(40):
        rng = random.Random(475101+seed)
        controls = []
        for _ in range(4):
            data = [rng.uniform(-1,1) for _ in range(4)]
            magnitude = math.hypot(*data)
            controls.append([x/magnitude for x in data])
        for scale in (1.,2.):
            # Nonunit accurate SLERP observations are compatibility only;
            # the supported unit domain is not expanded.
            data = [[scale*x for x in c] for c in controls]
            a = [previous.Quaternion(c[:]) for c in data]
            b = [q.Quaternion(c[:]) for c in data]
            for t in (-.25,.125,.37,.5,.875,1.25):
                for old,new in ((previous.quat_slerp(*a[:2],t),q.quat_slerp(*b[:2],t)),
                                (previous.squad4(*a,t),q.squad4(*b,t))):
                    assert packed(old.data) == packed(new.data)
                    interior += 1
                assert packed(previous.quat_squad(*a[:3],t).data) == packed(q.quat_squad(*b[:3],t).data)
                legacy += 1
            for t in (0.,1.):
                for old,new in ((previous.quat_slerp(*a[:2],t),q.quat_slerp(*b[:2],t)),
                                (previous.squad4(*a,t),q.squad4(*b,t))):
                    endpoints += 1
                    changed_endpoints += packed(old.data) != packed(new.data)
                    largest_endpoint_difference = max(largest_endpoint_difference,
                                                      *(abs(x-y) for x,y in zip(old.data,new.data)))
    characterized = []
    for data in ([float('nan'),0.,1.,0.],[float('inf'),0.,1.,0.],
                 [-float('inf'),0.,1.,0.],[1.,0.,0.,0.],[1.,2.,3.]):
        for t in (0.,1.,.37,float('nan'),float('inf')):
            outcomes = []
            for module in (previous,q):
                try:
                    values = module.quat_slerp(module.Quaternion(list(data)),module.Quaternion(),t).data
                    outcomes.append([('nan' if math.isnan(x) else x) for x in values])
                except (TypeError,IndexError,ValueError,ZeroDivisionError) as error:
                    outcomes.append({'exception':type(error).__name__, 'message':str(error)})
            assert outcomes[0] == outcomes[1], (data,t,outcomes)
            characterized.append({'input':data,'t':t,'outcome':outcomes[1]})
    return {'interior_slerp_squad4_bit_identical_calls':interior,
            'endpoint_calls':endpoints, 'changed_endpoint_calls':changed_endpoints,
            'maximum_endpoint_absolute_difference':largest_endpoint_difference,
            'legacy_squad_bit_identical_calls':legacy,
            'nonunit_observations_not_a_supported_domain':True,
            'invalid_nonfinite_characterizations':characterized}


def summary(values):
    median = statistics.median(values)
    return {'median':median,'mad':statistics.median(abs(x-median) for x in values),
            'min':min(values),'max':max(values)}


def axis(angle):
    return [math.cos(angle),0.,0.,math.sin(angle)]


def json_finite(value):
    if isinstance(value,float) and not math.isfinite(value):
        return 'nan' if math.isnan(value) else 'infinity' if value > 0 else '-infinity'
    if isinstance(value,dict): return {k:json_finite(v) for k,v in value.items()}
    if isinstance(value,list): return [json_finite(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--trials',type=int,default=9)
    parser.add_argument('--iterations',type=int,default=20000)
    args = parser.parse_args()
    if args.trials < 3 or args.iterations < 100:
        parser.error('Use at least three trials and 100 iterations')
    previous,old_source = old_module()
    comparisons = compatibility(previous)
    unit = math.ulp(0.)
    definitions = [
        ('slerp_raw_ordinary','raw', [axis(.1),axis(.9)],.37),
        ('slerp_method_ordinary','method',[axis(.1),axis(.9)],.37),
        ('slerp_method_close','method',[axis(.1),axis(.100001)],.37),
        ('slerp_method_tiny_normal','method',[axis(0.),[1.,1e-200,0.,0.]],.37),
        ('slerp_method_identical','method',[axis(.4),axis(.4)],.37),
        ('slerp_method_antipodal','method',[axis(.4),[-x for x in axis(.4)]],.37),
        ('slerp_endpoint_0','method',[axis(.1),axis(.9)],0.),
        ('slerp_endpoint_1','method',[axis(.1),axis(.9)],1.),
        ('slerp_subnormal_endpoint','method',[[1.,0.,0.,0.],[1.,unit,0.,0.]],1.),
        ('slerp_subnormal_interior','method',[[1.,0.,0.,0.],[1.,3*unit,0.,0.]],.375),
        ('squad4_ordinary','squad4',[axis(.1),axis(.9),axis(.3),axis(1.1)],.37),
        ('squad4_endpoint_0','squad4',[axis(.1),axis(.9),axis(.3),axis(1.1)],0.),
        ('squad4_endpoint_1','squad4',[axis(.1),axis(.9),axis(.3),axis(1.1)],1.),
        ('legacy_squad_control','legacy',[axis(.1),axis(.9),axis(.3)],.37),
    ]
    workloads = []
    for name,operation,data,t in definitions:
        work = {'name':name,'operation':operation,'inputs':data,'t':t,
                'operations_per_trial':args.iterations,'before':[],'after':[],'paired_ratios':[]}
        for side,module in (('before',previous),('after',q)):
            objects = [module.Quaternion(values[:]) for values in data]
            if operation == 'raw':
                work[side+'_call'] = lambda f=module.quat_slerp,controls=objects,t=t:f(*controls,t)
            elif operation == 'method':
                work[side+'_call'] = lambda controls=objects,t=t:controls[0].slerp(controls[1],t)
            else:
                function = module.squad4 if operation == 'squad4' else module.quat_squad
                work[side+'_call'] = lambda f=function,controls=objects,t=t:f(*controls,t)
            for _ in range(200): work[side+'_call']()
        workloads.append(work)
    for trial in range(args.trials):
        order = list(workloads)
        random.Random(475102+trial).shuffle(order)
        for work in order:
            pair = {}
            for side in (('before','after') if trial%2 == 0 else ('after','before')):
                was_enabled = gc.isenabled(); gc.disable()
                try:
                    call = work[side+'_call']
                    start = time.process_time_ns()
                    for _ in range(args.iterations): call()
                    elapsed = (time.process_time_ns()-start)/args.iterations
                finally:
                    if was_enabled: gc.enable()
                work[side].append(elapsed); pair[side] = elapsed
            work['paired_ratios'].append(pair['after']/pair['before'])
    for work in workloads:
        work['before_ns'] = summary(work['before']); work['after_ns'] = summary(work['after'])
        work['paired_after_before'] = summary(work['paired_ratios'])
        del work['before_call'],work['after_call']
    result = {'schema':1,'base':BASE,'compatibility':comparisons,
              'source_sha256':{'before':hashlib.sha256(old_source.encode()).hexdigest(),
                               'after':hashlib.sha256((ROOT/'gem/quaternion.py').read_bytes()).hexdigest()},
              'environment':{'python':sys.version,'platform':platform.platform(),'machine':platform.machine()},
              'method':{'clock':'process_time_ns','trials':args.trials,'iterations':args.iterations,
                        'warmup_calls_per_side':200,'setup_excluded':True,'gc_disabled_only_during_timing':True,
                        'loop_and_allocation_cost_included':True,
                        'order':'alternating implementations; deterministic shuffle seed 475102+trial'},
              'workloads':workloads}
    args.output.write_text(json.dumps(json_finite(result),indent=2,allow_nan=False)+'\n')
    for work in workloads:
        print('{}: {:.3f} -> {:.3f} us; paired ratio {:.3f}'.format(work['name'],
              work['before_ns']['median']/1000,work['after_ns']['median']/1000,
              work['paired_after_before']['median']))


if __name__ == '__main__':
    main()
