"""Minimal geometry reproductions; runs against source or an isolated install.

PYTHONPATH=. python audit/geometry_reproductions.py --output /tmp/geometry.json
venv/bin/python -I audit/geometry_reproductions.py --installed --output /tmp/wheel.json
"""
import argparse
import hashlib
import importlib.metadata
import json
import math
from pathlib import Path
import platform
import sys

from gem import matrix, plane, quaternion, ray, vector


def V(values):
    return vector.Vector(len(values),list(values))


def clean(value):
    if isinstance(value,float) and not math.isfinite(value):return str(value)
    if isinstance(value,list):return [clean(v) for v in value]
    return value


def capture(function):
    try:return {'value':clean(function())}
    except (ZeroDivisionError,OverflowError,ValueError) as error:
        return {'exception':type(error).__name__,'message':str(error)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--installed',action='store_true')
    args=parser.parse_args()
    root=Path(__file__).resolve().parents[1]
    modules=[vector,plane,ray,matrix,quaternion]
    if args.installed:
        assert all(Path(m.__file__).is_relative_to(sys.prefix) for m in modules)
        assert all(not Path(m.__file__).is_relative_to(root) for m in modules)
    # Independent literal geometry and unit basis answers, before findings.
    p=plane.Plane();p.fromCoeffs(2,3,6,-26)
    assert p.dot(V([1,2,3,1]))==0
    n=p.normalize();assert n.normal.vector==[2/7,3/7,6/7]
    assert n.d==-26/7
    r=ray.Ray(V([1,2,3]),V([0,0,5]));copy=r.duplicate()
    assert copy.distance==5 and copy.start.vector is not r.start.vector
    q=quaternion.Quaternion([.5,.5,.5,.5]);qstore=q.data
    assert r.rotateUsingQuaternion(q) is None
    assert r.start.vector==[3.,1.,2.] and r.dir.vector==[1.,0.,0.]
    assert q.data is qstore and q.data==[.5,.5,.5,.5]
    m=matrix.Matrix(4,[[1,0,0,0],[0,1,0,0],[0,0,1,0],[2,-3,4,1]])
    r.translate(m);assert r.start.vector==[5.,-2.,6.] and r.dir.vector==[1.,0.,0.]
    assert r.distance==5 and r.end.vector==[0.,0.,0.]
    findings=[]
    incident=[1.,-1e-9,0.];normal=[0.,1.,0.]
    actual=capture(lambda:vector.refract(1.,V(incident),V(normal)).vector)
    findings.append({'id':'4G4-A01','source':'gem/vector.py:65-71',
                     'inputs':{'IOR':1.,'incident':incident,'normal':normal},
                     'expected':incident,'actual':actual,
                     'reproduced':actual != {'value':incident}})
    o=float(2**52)
    for oblique in [False,True]:
        points=[[o,o,o],[o+1,o,o+float(oblique)],[o,o+1,o+float(oblique)]]
        expected=[-1/math.sqrt(3),-1/math.sqrt(3),1/math.sqrt(3)] if oblique else [0.,0.,1.]
        actual=capture(lambda:plane.Plane().bestFitNormal(list(map(V,points))).vector)
        findings.append({'id':'4G4-A02','source':'gem/plane.py:105-109',
                         'inputs':{'points':points},'exact_area_vector':[-int(oblique),-int(oblique),1],
                         'expected':expected,'actual':actual,
                         'reproduced':actual != {'value':expected}})
    observations=[
        {'classification':'documented limitation','name':'unscaled barycentric products',
         'expected':[.5,.25,.25],
         'actual':capture(lambda:V([.25e-100,.25e-100]).barycentric(V([0,0]),V([1e-100,0]),V([0,1e-100])))},
        {'classification':'documented limitation','name':'unscaled three-point cross underflow',
         'expected':[0.,0.,1.,0.],
         'actual':capture(lambda:plane.Plane().fromPoints(V([0,0,0]),V([1e-200,0,0]),V([0,1e-200,0])))},
        {'classification':'documented limitation','name':'plane normal magnitude outside binary64 range',
         'expected':[1/math.sqrt(3)]*4,
         'actual':capture(lambda:list(plane.normalize([1e308]*4)))},
        {'classification':'documented arithmetic limitation','name':'unscaled reflection product overflow',
         'expected':[-1e308,0.,0.],
         'actual':capture(lambda:vector.reflect(V([1e308,0,0]),V([1,0,0])).vector)},
        {'classification':'accuracy-policy question','name':'extreme polygon mean accumulation',
         'expected':1e308,
         'actual':capture(lambda:plane.Plane().bestFitD([V([0,0,1e308]),V([1,0,1e308])],V([0,0,1])))},
        {'classification':'accuracy-policy question','name':'enormous IOR at normal incidence',
         'expected':[0.,-1.,0.],
         'actual':capture(lambda:vector.refract(float(2**54),V([0,-1,0]),V([0,1,0])).vector)},
    ]
    dependencies={}
    for name in ['six','pytest']:
        try:dependencies[name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:pass
    result={'schema':'gem-geometry-reproductions-v1','environment':{'python':sys.version,
            'implementation':platform.python_implementation(),'platform':platform.platform(),
            'architecture':platform.machine(),'dependencies':dependencies},
            'installed_check':args.installed,'ordinary_smoke_passed':True,
            'runtime_modules':{m.__name__:{'path':m.__file__,'sha256':hashlib.sha256(Path(m.__file__).read_bytes()).hexdigest()} for m in modules},
            'findings':findings,'observations':observations,
            'limitations':'Observations outside established stability guarantees do not become new strict contract tests.'}
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'ordinary_smoke_passed':True,'installed_check':args.installed,
                      'findings_reproduced':[f['id'] for f in findings if f['reproduced']]}))


if __name__=='__main__':main()
