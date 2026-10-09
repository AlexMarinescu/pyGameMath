"""Paired Phase 4G-2R timings and alternative-evaluation evidence.

Run from a checkout; git supplies the immutable pre-repair implementation.
Only standard-library modules and the supported gem core are required.
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
sys.path.insert(0, str(ROOT))
from gem import bezier, legendre, spherical_harmonics as sh
from gem.vector import Vector

BASE = 'a52c8559ae9bcac6f86d84a9a6ef98b1b2274ac4'


def previous_module(base, name):
    source = subprocess.check_output(['git', 'show', base+':gem/'+name+'.py'], cwd=ROOT, text=True)
    module = types.ModuleType(name+'_before_repairs')
    exec(compile(source, base+':gem/'+name+'.py', 'exec'), module.__dict__)
    return module


def packed(values):
    return [struct.pack('!d', x) for x in values]


def casteljau(t, values):
    values = list(values)
    while len(values) > 1:
        values = [(1-t)*a+t*b for a, b in zip(values, values[1:])]
    return values[0]


def difference_casteljau(t, values):
    values = list(values)
    while len(values) > 1:
        values = [a+t*(b-a) for a, b in zip(values, values[1:])]
    return values[0]


def ordered_ordinary(degree, x):
    previous, current = 1., x
    if not degree:
        return previous
    for index in range(2, degree+1):
        previous, current = current, x*((2*index-1)/index)*current-((index-1)/index)*previous
    return current


def compatibility(before_b, before_l):
    counts = {'bezier_bitwise': 0, 'legendre_bitwise': 0, 'sh_bitwise': 0,
              'weighted_casteljau_changed_results': 0, 'ordered_legendre_changed_results': 0}
    for seed in range(16):
        rng = random.Random(472100+seed)
        for degree in [2, 3]:
            old = before_b.quadraticBezierPoint if degree == 2 else before_b.cubicBezierPoint
            new = bezier.quadraticBezierPoint if degree == 2 else bezier.cubicBezierPoint
            for dimension in [None, 2, 3, 4, 8]:
                data = [[rng.uniform(-4,4) for _ in range(dimension or 1)] for _ in range(degree+1)]
                controls = [row[0] for row in data] if dimension is None else [Vector(dimension,row) for row in data]
                for t in [-.5, 0., .2, .375, .5, .9, 1., 1.5]:
                    a, b = old(t,*controls), new(t,*controls)
                    a, b = ([a],[b]) if dimension is None else (a.vector,b.vector)
                    assert packed(a) == packed(b)
                    counts['bezier_bitwise'] += 1
                    candidate = [casteljau(t,[row[j] for row in data]) for j in range(dimension or 1)]
                    if packed(candidate) != packed(a):
                        counts['weighted_casteljau_changed_results'] += 1
    for degree in range(33):
        for order in sorted({0,min(1,degree),degree//2,degree}):
            for x in [-1.,-.875,-.37,0.,.2,.75,1.] + ([-2.,-1.25,1.25,2.] if order == 0 else []):
                a, b = before_l.Legendre(degree,order,x).run(), legendre.Legendre(degree,order,x).run()
                assert packed([a]) == packed([b])
                counts['legendre_bitwise'] += 1
                if order == 0 and packed([ordered_ordinary(degree,x)]) != packed([a]):
                    counts['ordered_legendre_changed_results'] += 1
    original = sh.Legendre
    try:
        for bands in [1,3,7,13]:
            for theta, phi in [(.2,-.3),(.73,1.27),(2.2,3.1),(0.,0.),(math.pi,.9)]:
                sh.Legendre = before_l.Legendre
                before = [sh.SPH(l,m,theta,phi) for l in range(bands) for m in range(-l,l+1)]
                sh.Legendre = original
                after = [sh.SPH(l,m,theta,phi) for l in range(bands) for m in range(-l,l+1)]
                assert packed(before) == packed(after)
                counts['sh_bitwise'] += len(after)
    finally:
        sh.Legendre = original
    return counts


def summarize(values):
    median = statistics.median(values)
    return {'median':median, 'mad':statistics.median(abs(x-median) for x in values),
            'min':min(values), 'max':max(values)}


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
    before_b, before_l = previous_module(args.base,'bezier'), previous_module(args.base,'legendre')
    checks = compatibility(before_b,before_l)
    workloads = []
    for degree in [2,3]:
        old = before_b.quadraticBezierPoint if degree == 2 else before_b.cubicBezierPoint
        new = bezier.quadraticBezierPoint if degree == 2 else bezier.cubicBezierPoint
        for dimension in [None,2,3,4]:
            for kind in ['ordinary','tiny']:
                t = .37 if kind == 'ordinary' else (1e-200 if degree == 2 else 1e-150)
                data = [(-1.)**i*(i+1) for i in range(degree+1)] if kind == 'ordinary' else [0.]*degree+[1e300]
                controls = data if dimension is None else [Vector(dimension,[p*(-1)**j for j in range(dimension)]) for p in data]
                calls = {'before':lambda f=old,p=controls,t=t:f(t,*p),
                         'after':lambda f=new,p=controls,t=t:f(t,*p)}
                if dimension is None:
                    calls['casteljau'] = lambda p=data,t=t:casteljau(t,p)
                else:
                    calls['casteljau'] = lambda p=controls,t=t,n=dimension:Vector(n,[casteljau(t,[v.vector[j] for v in p]) for j in range(n)])
                workloads.append({'name':'bezier_{}_{}_{}'.format(degree,dimension or 'scalar',kind),
                                  'input':{'t':t,'scalar_controls':data,'dimension':dimension},'calls':calls})
    for degree,order,x in [(0,0,.37),(2,0,.37),(2,2,.37),(12,0,.37),(12,5,.37),
                            (32,12,.875),(2,0,2.),(12,0,2.),(2,0,1e154),(3,0,4e102)]:
        name='legendre_{}_{}_{}'.format(degree,order,x)
        workloads.append({'name':name,'input':{'degree':degree,'order':order,'x':x},
                          'calls':{'before':lambda l=degree,m=order,x=x:before_l.Legendre(l,m,x).run(),
                                   'after':lambda l=degree,m=order,x=x:legendre.Legendre(l,m,x).run()}})
        if degree == 2 and order == 0:
            old_object,new_object = before_l.Legendre(degree,order,x),legendre.Legendre(degree,order,x)
            workloads.append({'name':name+'_bound_run','input':{'degree':degree,'order':order,'x':x},
                              'calls':{'before':old_object.run,'after':new_object.run}})
    for work in workloads:
        work['samples_ns']={side:[] for side in work['calls']}
        work['paired_after_before']=[]
    enabled=gc.isenabled();gc.disable()
    try:
        for work in workloads:
            for call in work['calls'].values():
                for _ in range(200):call()
        for trial in range(args.trials):
            order=list(workloads);random.Random(472101+trial).shuffle(order)
            for work in order:
                sides=list(work['calls']) if trial%2==0 else list(reversed(work['calls']))
                pair={}
                for side in sides:
                    pair[side]=latency(work['calls'][side],args.iterations)
                    work['samples_ns'][side].append(pair[side])
                work['paired_after_before'].append(pair['after']/pair['before'])
    finally:
        if enabled:gc.enable()
    for work in workloads:
        del work['calls']
        work['summary_ns']={side:summarize(values) for side,values in work['samples_ns'].items()}
        work['paired_summary']=summarize(work['paired_after_before'])
        print('{}: {:.3f} -> {:.3f} us; ratio {:.3f}'.format(work['name'],
              work['summary_ns']['before']['median']/1000,work['summary_ns']['after']['median']/1000,
              work['paired_summary']['median']))
    alternative_evidence=[]
    for t,values in [(1e-200,[0.,0.,1e300]),(1e-150,[0.,0.,0.,1e300]),(.5,[-1e308,1e308,-1e308])]:
        actual=bezier.quadraticBezierPoint(t,*values) if len(values)==3 else bezier.cubicBezierPoint(t,*values)
        alternative_evidence.append({'t':t,'controls':values,'selected':actual,
                                     'weighted_casteljau':casteljau(t,values),
                                     'difference_casteljau':str(difference_casteljau(t,values))})
    output={'schema':1,'baseline_commit':args.base,'compatibility':checks,
            'source_sha256':{name:hashlib.sha256((ROOT/'gem'/name).read_bytes()).hexdigest() for name in ['bezier.py','legendre.py']},
            'environment':{'python':sys.version,'implementation':platform.python_implementation(),
                           'platform':platform.platform(),'machine':platform.machine()},
            'methodology':{'clock':'process_time_ns','trials':args.trials,'operations_per_trial':args.iterations,
                           'setup_excluded':True,'warmup_calls_per_side':200,'gc_disabled':True,
                           'loop_overhead_subtracted':False,'shuffle_seed':'472101+trial','order':'alternating'},
            'workloads':workloads,'alternative_evidence':alternative_evidence}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(output,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
