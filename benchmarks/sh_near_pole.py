"""Matched SH repair timings, ordinary rounding, and HDR reference differences.

Standard-library development tool; immutable pre-repair source comes from git.
Run after regenerating HDR outputs into a separate directory, never the goldens.
"""
import argparse
import cProfile
import gc
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import re
import statistics
import struct
import subprocess
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from gem import spherical_harmonics as sh
from gem.quaternion import Quaternion
from gem.vector import Vector
from benchmarks.verify_hdr_reference import png_rgb

BASE = '72fb9dac9beadc6c7ca1eae8a615b4ad86441254'


def before_module(base):
    source = subprocess.check_output(['git','show',base+':gem/spherical_harmonics.py'],cwd=ROOT,text=True)
    module = types.ModuleType('sh_before_near_pole_repair')
    exec(compile(source,base+':gem/spherical_harmonics.py','exec'),module.__dict__)
    return module,source


def rows(count=9):
    return [[(-1)**i*(i+.25)/8,(i+.5)/16,0.] for i in range(count)]


def compared(a,b):
    values_a = [a] if isinstance(a,(int,float)) else a
    values_b = [b] if isinstance(b,(int,float)) else b
    return [(struct.pack('!d',x)==struct.pack('!d',y),abs(x-y)) for x,y in zip(values_a,values_b)]


def compatibility(before):
    metrics = {}
    def record(name,a,b):
        m = metrics.setdefault(name,{'components':0,'bit_identical':0,'maximum_absolute_difference':0.})
        for same,delta in compared(a,b):
            m['components'] += 1; m['bit_identical'] += same
            m['maximum_absolute_difference'] = max(m['maximum_absolute_difference'],delta)
            assert delta <= 2e-13*max(1.,abs(a) if isinstance(a,(int,float)) else max(map(abs,a)))
    for theta in [.1,.3,.73,1.2,math.pi/2,2.1,math.pi-.3]:
        for phi in [-1.,0.,.7,1.27]:
            for l in range(13):
                for m in range(-l,l+1):
                    record('ordinary_SPH',before.SPH(l,m,theta,phi),sh.SPH(l,m,theta,phi))
            for bands in [1,3,5,9]:
                record('ordinary_basis',before._basis(bands,theta,phi),sh._basis(bands,theta,phi))
    rng = random.Random(473101)
    for _ in range(128):
        d = [rng.uniform(-1,1) for _ in range(3)]; n = math.hypot(*d); d = [v/n for v in d]
        for count in [1,4,9,25]:
            c = [[rng.uniform(-2,2) for _ in range(3)] for _ in range(count)]
            record('ordinary_reconstruct',before.reconstruct(c,d),sh.reconstruct(c,d))
    for theta in [0.,math.pi]:
        for l in range(32):
            for m in range(-l,l+1):
                a,b=before.SPH(l,m,theta,.73),sh.SPH(l,m,theta,.73)
                assert struct.pack('!d',a)==struct.pack('!d',b)
                record('exact_poles',a,b)
    c = rows(); q = Quaternion([math.sqrt(.5),0.,0.,math.sqrt(.5)])
    for function in ['rotate_coefficients','convolve_diffuse','legacy_to_canonical']:
        arguments = (c,q) if function == 'rotate_coefficients' else (c,)
        a,b = getattr(before,function)(*arguments),getattr(sh,function)(*arguments)
        assert a == b
        for x,y in zip(a,b): record(function,x,y)
    # Invalid and out-of-domain historical cases are observations, not new APIs.
    domains = []
    for l,m,theta,phi in [(2,1,-.7,.3),(2,1,-1e-9,.3),(2,1,2*math.pi,.3),
                         (2,1,math.pi+.1,.3),(2,1,math.nan,.3),(2,1,math.inf,.3),
                         (2,1,None,.3),(2,1.5,.7,.3),(2,3,.7,.3),(0,0,math.nan,.3)]:
        def capture(module):
            try: return {'value_hex':module.SPH(l,m,theta,phi).hex()}
            except Exception as e: return {'exception':type(e).__name__}
        old,new=capture(before),capture(sh); assert old == new
        domains.append({'arguments':repr((l,m,theta,phi)),'before':old,'after':new})
    return {'metrics':metrics,'historical_domain_observations':domains,
            'limitations':'Ordinary last-bit differences are recorded; independent mathematical tests remain the correctness oracle.'}


def json_differences(a,b,path=''):
    if isinstance(a,dict):
        assert isinstance(b,dict) and a.keys() == b.keys()
        return [d for k in a for d in json_differences(a[k],b[k],path+'/'+k)]
    if isinstance(a,list):
        assert isinstance(b,list) and len(a) == len(b)
        return [d for i,(x,y) in enumerate(zip(a,b)) for d in json_differences(x,y,path+'/'+str(i))]
    if isinstance(a,(int,float)) and isinstance(b,(int,float)):
        if a == b: return []
        delta=abs(a-b)
        assert delta <= 3e-14*max(1.,abs(a),abs(b)),path
        return [{'path':path,'reference':a,'generated':b,'absolute_difference':delta}]
    assert a == b,(path,a,b)
    return []


def hdr_comparison(directory):
    assert directory.resolve() != (ROOT/'examples/output').resolve()
    files = {}
    for reference in sorted((ROOT/'examples/output').iterdir()):
        generated = directory/reference.name
        old,new = reference.read_bytes(),generated.read_bytes()
        record = {'reference_sha256':hashlib.sha256(old).hexdigest(),
                  'generated_sha256':hashlib.sha256(new).hexdigest(),'bytes_identical':old==new}
        if reference.suffix == '.json':
            changes = json_differences(json.loads(old),json.loads(new))
            record.update(numeric_changes=changes,changed_numeric_fields=len(changes),
                          maximum_absolute_difference=max((d['absolute_difference'] for d in changes),default=0.))
        elif reference.suffix == '.png':
            ow,oh,op=png_rgb(reference);nw,nh,np=png_rgb(generated)
            assert (ow,oh)==(nw,nh)
            record.update(width=nw,height=nh,decoded_pixels_identical=op==np,
                          changed_channel_bytes=sum(a!=b for a,b in zip(op,np)),
                          maximum_channel_byte_difference=max(abs(a-b) for a,b in zip(op,np)))
        elif reference.suffix == '.glsl':
            pattern=r'vec3\(([^)]*)\)'
            floats=lambda data:[float(v.strip()) for row in re.findall(pattern,data.decode()) for v in row.split(',')]
            a,b=floats(old),floats(new);assert len(a)==len(b)==27
            differences=[abs(x-y) for x,y in zip(a,b)]
            assert max(differences)<=3e-14*max(1.,max(map(abs,a)))
            record.update(changed_numeric_constants=sum(x!=y for x,y in zip(a,b)),maximum_absolute_difference=max(differences))
        files[reference.name] = record
    return {'regeneration_command':f'python -m examples.hdr_sh.regenerate --output-dir {directory}',
            'files':files,'goldens_modified':False,
            'limitations':'Same-host byte and numerical comparisons; floating/libm/compression bytes are not a cross-platform correctness guarantee.'}


def workloads(module,iterations):
    c = rows(); c25 = rows(25); q = Quaternion([math.sqrt(.5),0,0,math.sqrt(.5)])
    samples=[]
    for i in range(16):
        theta=.1+(math.pi-.2)*(i+.5)/16;phi=.2+i*.7
        s=module.SPHSample(theta,phi,Vector(3,[0,0,1]),9);s.values=sh._basis(3,theta,phi);samples.append(s)
    colors=[[1.,2.,.5] for _ in samples]
    image=[[[1.+j/32,2.+i/16,.5] for j in range(32)] for i in range(16)]
    cases = {
        'SPH_L0':lambda:module.SPH(0,0,.73,1.27),
        'SPH_L1_ordinary':lambda:module.SPH(1,1,.73,1.27),
        'SPH_L1_near_north':lambda:module.SPH(1,1,1e-9,.3),
        'SPH_L2_ordinary':lambda:module.SPH(2,-1,.73,1.27),
        'SPH_L2_near_south':lambda:module.SPH(2,-1,math.pi-1e-9,.3),
        'SPH_L12_m6':lambda:module.SPH(12,6,.73,1.27),
        'basis_B3_ordinary':lambda:module._basis(3,.73,1.27),
        'basis_B3_near':lambda:module._basis(3,1e-9,.3),
        'basis_B9_ordinary':lambda:module._basis(9,.73,1.27),
        'reconstruct_L2_ordinary':lambda:module.reconstruct(c,[.3,.4,math.sqrt(.75)]),
        'reconstruct_L2_near_north':lambda:module.reconstruct(c,[1e-9,-2e-9,1.]),
        'reconstruct_L2_near_south':lambda:module.reconstruct(c,[1e-100,-2e-100,-1.]),
        'reconstruct_B5':lambda:module.reconstruct(c25,[.3,.4,math.sqrt(.75)]),
        'rotate_L2_RGB_control':lambda:module.rotate_coefficients(c,q),
        'convolve_L2_control':lambda:module.convolve_diffuse(c),
        'project_samples_16_control':lambda:module.project_radiance(samples,colors),
        'project_angular_32x16':lambda:module.project_angular_probe(image),
        'GenerateSamples_4x4_B3':lambda:module.GenerateSamples(4,3),
    }
    return {name:(function,max(1,iterations//(512 if name.startswith('project_angular') else
                        32 if name.startswith('GenerateSamples') else 4 if name.startswith('basis_B9') else 1)))
            for name,function in cases.items()}


def timed(function,count):
    active=gc.isenabled();gc.disable()
    try:
        start=time.process_time_ns()
        for _ in range(count): function()
        return (time.process_time_ns()-start)/count
    finally:
        if active:gc.enable()


def summary(values):
    median=statistics.median(values)
    return {'median_ns':median,'mad_ns':statistics.median(abs(v-median) for v in values),
            'minimum_ns':min(values),'maximum_ns':max(values),'samples_ns':values}


def profile(function):
    p=cProfile.Profile();p.enable()
    for _ in range(1000):function()
    p.disable()
    return sorted([{'function':entry.code if isinstance(entry.code,str) else entry.code.co_name,
                    'calls':entry.callcount,'total_seconds':entry.totaltime,'self_seconds':entry.inlinetime}
                   for entry in p.getstats()],key=lambda row:row['self_seconds'],reverse=True)[:10]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base',default=BASE);parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--hdr-directory',type=Path,required=True)
    parser.add_argument('--trials',type=int,default=9);parser.add_argument('--iterations',type=int,default=20000)
    args=parser.parse_args();before,source=before_module(args.base)
    correctness=compatibility(before);hdr=hdr_comparison(args.hdr_directory)
    old,new=workloads(before,args.iterations),workloads(sh,args.iterations)
    assert old.keys()==new.keys()
    data={name:{'before':[],'after':[]} for name in old}
    rng_state=random.getstate()
    try:
        random.seed(473103)
        for cases in [old,new]:
            for function,_ in cases.values():
                for _ in range(100):function()
        for trial in range(args.trials):
            names=list(old);random.Random(473102+trial).shuffle(names)
            for name in names:
                order=[('before',old),('after',new)] if trial%2==0 else [('after',new),('before',old)]
                for label,cases in order:
                    function,count=cases[name]
                    if name.startswith('GenerateSamples'): random.seed(473104+trial)
                    data[name][label].append(timed(function,count))
    finally: random.setstate(rng_state)
    result={}
    for name,values in data.items():
        ratios=[b/a for a,b in zip(values['before'],values['after'])]
        result[name]={'calls_per_trial':old[name][1],'before':summary(values['before']),
                      'after':summary(values['after']),'paired_after_before':summary(ratios)}
    # Profiling is separate from latency measurement.
    profiles={label:profile(cases['basis_B3_ordinary'][0]) for label,cases in [('before',old),('after',new)]}
    cpu=next((line.split(':',1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')),None)
    output={'schema':'gem-sh-near-pole-repair-v1','base':args.base,
            'source_sha256':{'before':hashlib.sha256(source.encode()).hexdigest(),
                             'after':hashlib.sha256((ROOT/'gem/spherical_harmonics.py').read_bytes()).hexdigest()},
            'environment':{'python':sys.version,'implementation':platform.python_implementation(),
                           'platform':platform.platform(),'architecture':platform.machine(),'cpu':cpu},
            'methodology':{'trials':args.trials,'iterations':args.iterations,'clock':'process_time_ns',
                           'warmup_calls_per_workload_and_implementation':100,'gc_during_timing':False,
                           'trial_order':'alternating before/after, shuffled workloads seed 473102+trial',
                           'setup_excluded':True,'loop_and_return_allocation_included':True,
                           'random_generation':'Warmup seed 473103; each paired generation trial resets seed 473104+trial outside timing; restored on exit; comparison seed 473101',
                           'latency_ratio':'median of paired after/before; not the ratio of independent medians',
                           'limitations':'Shared-host frequency/load not controlled; medians/MADs describe this run, not universal guarantees.'},
            'compatibility':correctness,'hdr_reference':hdr,'workloads':result,'profiles':profiles}
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')
    for name,row in result.items():
        print(f'{name}: {row["before"]["median_ns"]/1000:.3f} -> {row["after"]["median_ns"]/1000:.3f} us; paired after/before {row["paired_after_before"]["median_ns"]:.3f}')


if __name__=='__main__':main()
