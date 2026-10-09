"""Bit-preservation and unchanged-code checks against merged Phase 3C."""
import ast
import hashlib
import json
import math
from pathlib import Path
import random
import struct
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from benchmarks.bezier_sh import baseline, cases, BASE
from gem import bezier, spherical_harmonics as sh
from gem.vector import Vector
from gem.quaternion import Quaternion


def bits(value):
    if isinstance(value,Vector):return (value.size,bits(value.vector))
    if isinstance(value,(tuple,list)):return tuple(bits(v) for v in value)
    return struct.pack('!d',float(value))


def main():
    (ob,os),hashes=baseline();counts={};rng=random.Random(304)
    def same(name,a,b):
        assert bits(a)==bits(b),name
        counts[name]=counts.get(name,0)+1
    for dimension in (0,1,2,3,4,8):
        for _ in range(50):
            coords=[[rng.uniform(-10,10) for _ in range(dimension)] for _ in range(4)]
            controls=[Vector(dimension,p) for p in coords];t=rng.uniform(-.5,1.5)
            for name,degree in (('quadraticBezierPoint',2),('cubicBezierPoint',3)):
                same(name,getattr(ob,name)(t,*controls[:degree+1]),getattr(bezier,name)(t,*controls[:degree+1]))
    for dimension in (1,2,3):
        for _ in range(50):
            controls=[tuple(rng.uniform(-4,4) for _ in range(dimension)) for _ in range(4)]
            same('flatness',max(ob._chord_distance(p,controls[0],controls[-1]) for p in controls[1:-1]),bezier._flatness(controls))
            same('subdivide',ob._subdivide(controls,.01),bezier._subdivide(controls,.01))
    for bands in (1,2,3,5,9,13):
        for _ in range(30):
            theta=rng.uniform(0,math.pi);phi=rng.uniform(-math.pi,math.pi)
            same('basis',os._basis(bands,theta,phi),sh._basis(bands,theta,phi))
    for channels in (1,3):
        for _ in range(100):
            rows=[[rng.uniform(-4,4) for _ in range(channels)] for _ in range(9)]
            coefficients=[r[0] for r in rows] if channels==1 else rows
            q=Quaternion([rng.uniform(-1,1) for _ in range(4)]).normalize()
            same('rotation',os.rotate_coefficients(coefficients,q),sh.rotate_coefficients(coefficients,q))
    before,_=cases(ob,os);after,_=cases(bezier,sh)
    for name in before:same('benchmark_'+name,before[name](),after[name]())
    unchanged={}
    changed={'bezier':{'cubicBezierPoint','quadraticBezierPoint','_subdivide'},'spherical_harmonics':{'_basis','rotate_coefficients'}}
    added={'bezier':{'_flatness'},'spherical_harmonics':{'_basis_layout'}}
    for name in changed:
        old=subprocess.check_output(['git','show',BASE+':gem/'+name+'.py'],cwd=ROOT,text=True)
        new=(ROOT/'gem'/f'{name}.py').read_text()
        nodes=lambda source:{n.name:n for n in ast.parse(source).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        a,b=nodes(old),nodes(new)
        assert set(b)-set(a)==added[name]
        keys=[k for k in a if k not in changed[name]]
        for key in keys:assert ast.dump(a[key])==ast.dump(b[key]),(name,key)
        unchanged[name]=keys
    result={'baseline_commit':BASE,'bit_identical_checks':counts,'total':sum(counts.values()),
            'unchanged_definitions':unchanged,'baseline_sha256':hashes,
            'current_sha256':{n:hashlib.sha256((ROOT/'gem'/f'{n}.py').read_bytes()).hexdigest() for n in hashes}}
    Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
    print(f"{result['total']} bit-identical results; unrelated definitions unchanged")


if __name__=='__main__':main()
