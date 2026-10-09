"""Bit-preservation and unchanged-code checks against merged Phase 3B."""
import ast
import hashlib
import json
import math
from pathlib import Path
import random
import struct
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from benchmarks.vector_quaternion import baseline, BASE
from gem import vector, quaternion


def bits(values):
    return b''.join(struct.pack('!d',float(value)) for value in values)


def main():
    old,hashes=baseline();ov,oq=old
    counts={}
    def same(name,a,b):
        assert bits(a)==bits(b), (name,a,b)
        counts[name]=counts.get(name,0)+1
    rng=random.Random(303)
    for size in (0,1,2,3,4,8):
        for _ in range(100):
            data=[math.ldexp(rng.uniform(-1,1),rng.randint(-1074,1023)) for _ in range(size)]
            same('normalize',ov.normalize(size,data),vector.normalize(size,data))
    for data in ([0.,-0.,0.,-0.],[1e308]*4,[1e300,1e-300,-1.,0.],
                 [float('inf'),0.,1.,2.],[-float('inf'),1.,2.,3.],[float('nan'),0.,1.,2.]):
        same('normalize',ov.normalize(4,data),vector.normalize(4,data))
        same('quaternion_normalize',oq.quat_normalize(data),quaternion.quat_normalize(data))
        same('conjugate',oq.quat_conjugate(data),quaternion.quat_conjugate(data))
    for _ in range(200):
        data=[rng.uniform(-4,4) for _ in range(4)];point=[rng.uniform(-10,10) for _ in range(3)]
        same('conjugate',oq.quat_conjugate(data),quaternion.quat_conjugate(data))
        same('rotation',oq.quat_rotate_vector(oq.Quaternion(data),ov.Vector(3,point)).vector,
                        quaternion.quat_rotate_vector(quaternion.Quaternion(data),vector.Vector(3,point)).vector)
        controls=[[rng.uniform(-1,1) for _ in range(4)] for _ in range(4)]
        oldq=[oq.Quaternion(oq.quat_normalize(x)) for x in controls]
        newq=[quaternion.Quaternion(quaternion.quat_normalize(x)) for x in controls]
        t=rng.random()
        for name in ('quat_lerp','quat_slerp','quat_slerp_no_invert','quat_squad','squad4'):
            n=4 if name=='squad4' else 3 if name=='quat_squad' else 2
            same(name,getattr(oq,name)(*oldq[:n],t).data,getattr(quaternion,name)(*newq[:n],t).data)
    errors=[]
    invalid=[('normalize',(3,[1.,2.])),('normalize',(3,[10**1000,0,0])),('quat_conjugate',([1,2,3],))]
    for name,args in invalid:
        pair=[]
        for mod in ((ov,vector) if name=='normalize' else (oq,quaternion)):
            try:getattr(mod,name)(*args)
            except Exception as error:pair.append(type(error).__name__)
            else:raise AssertionError('Expected historical error')
        assert pair[0]==pair[1]
        errors.append({'operation':name,'exception':pair[0]})
    unchanged={}
    changed={'vector':{'normalize'},'quaternion':{'quat_conjugate','quat_rotate_vector','quat_lerp','quat_slerp','quat_slerp_no_invert'}}
    import subprocess
    for name in changed:
        old_text=subprocess.check_output(['git','show',BASE+':gem/'+name+'.py'],cwd=ROOT,text=True)
        new_text=(ROOT/'gem'/f'{name}.py').read_text()
        nodes=lambda text:{n.name:n for n in ast.parse(text).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
        before,after=nodes(old_text),nodes(new_text)
        assert set(after)-set(before)==({'_quat_blend'} if name=='quaternion' else set())
        names=[n for n in before if n not in changed[name]]
        for n in names:assert ast.dump(before[n])==ast.dump(after[n]),(name,n)
        unchanged[name]=names
    result={'baseline_commit':BASE,'bit_identical_checks':counts,'total':sum(counts.values()),
            'invalid_input_checks':errors,'unchanged_definitions':unchanged,'baseline_sha256':hashes,
            'current_sha256':{n:hashlib.sha256((ROOT/'gem'/f'{n}.py').read_bytes()).hexdigest() for n in hashes}}
    Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
    print(f"{result['total']} bit-identical outputs; unrelated definitions unchanged")


if __name__=='__main__':main()
