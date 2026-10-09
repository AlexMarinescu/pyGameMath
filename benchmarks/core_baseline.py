"""Deterministic, dependency-free benchmark/profiling harness (not gem runtime)."""
import argparse
import cProfile
import importlib.metadata
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import pstats
import statistics
import subprocess
import sys
import time
import timeit
import tracemalloc

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import bezier, matrix, plane, quaternion, ray, spherical_harmonics as sh, vector

HISTORICAL = '2cbd899a47546b8c40f8244f57799969c86b0d87'


def V(*values):
    return vector.Vector(len(values), list(values))


def path(tolerance, segments=1):
    p = bezier.BezierPath()
    points = [V(0, 0, 0)]
    for i in range(segments):
        points.extend([V(i + .1, 2, .5), V(i + .9, -2, -.5), V(i + 1, 0, 0)])
    p.setControlPoints(points)
    p.minimum_sqr_distance = tolerance * tolerance
    return p


def prepare():
    """All data generation and immutable inputs live outside timed calls."""
    calls, details = {}, {}
    def add(name, fn, domain='ordinary', **metadata):
        calls[name] = fn
        details[name] = {'domain': domain, **metadata}
    for n in (2, 3, 4):
        a = V(*[(-1)**i * (i + .25) for i in range(n)])
        b = V(*[i + .5 for i in range(n)])
        add(f'vector{n}_allocate', lambda n=n: vector.Vector(n))
        for label, fn in [('add', lambda a=a,b=b:a+b), ('subtract', lambda a=a,b=b:a-b),
                          ('scalar_multiply',lambda a=a:a*1.25),('dot',lambda a=a,b=b:a.dot(b)),
                          ('magnitude',a.magnitude),('normalize',a.normalize)]:
            add(f'vector{n}_{label}',fn)
        rows = [[3.0*(i==j) + .13*(i+1)/(j+1) for j in range(n)] for i in range(n)]
        m = matrix.Matrix(n, rows)
        other = matrix.Matrix(n, [[.25*(i-j)+(i==j) for j in range(n)] for i in range(n)])
        add(f'matrix{n}_allocate',lambda n=n:matrix.Matrix(n))
        add(f'matrix{n}_raw_multiply',lambda rows=rows:matrix.matrix_multiply(rows,rows))
        add(f'matrix{n}_multiply',lambda m=m,other=other:m*other)
        add(f'matrix{n}_vector',lambda m=m,a=a:m*a)
        add(f'matrix{n}_determinant',m.det)
        add(f'matrix{n}_inverse',m.inverse)
        add(f'matrix{n}_raw_inverse',lambda rows=rows,n=n:getattr(matrix,'inverse'+str(n))(rows))
        add(f'matrix{n}_transpose',m.transpose)
        add(f'vector{n}_transform',lambda a=a,rows=rows:vector.transform(a.size,a.vector,rows))
        if n>2:
            for scale in (1e-200,1e200):
                scaled=[[value*scale for value in row] for row in rows]
                add(f'matrix{n}_raw_inverse_scale_{scale:g}',
                    lambda scaled=scaled,n=n:getattr(matrix,'inverse'+str(n))(scaled),'extreme uniform scale')
    add('vector3_cross',lambda:vector.cross(V3,V3B))
    for values,label in [([1e300,-2e300,3e300],'large'),([1e-300,-2e-300,3e-300],'tiny'),
                         ([1e300,1e-300,-1],'mixed')]:
        v=V(*values)
        add('vector3_magnitude_'+label,v.magnitude,label)
        add('vector3_normalize_'+label,v.normalize,label)
    q=quaternion.quat_from_axis_angle(V(1,2,3),73)
    q1=quaternion.quat_from_axis_angle(V(-2,1,.5),-41)
    qnear=quaternion.quat_from_axis_angle(V(1,2,3),73.001)
    qm=q.toMatrix()
    for label,values in [('large',[1e300,-2e300,3e300,-4e300]),('tiny',[1e-300,-2e-300,3e-300,-4e-300])]:
        extreme=quaternion.Quaternion(values)
        add('quaternion_normalize_'+label,extreme.normalize,label)
    for label,fn in [('allocate',lambda:quaternion.Quaternion()),('multiply',lambda:q*q1),
                     ('normalize',q.normalize),('inverse',q.inverse),
                     ('rotate_vector',lambda:quaternion.quat_rotate_vector(q,V3)),
                     ('slerp',lambda:q.slerp(q1,.37)),('slerp_near',lambda:q.slerp(qnear,.37)),
                     ('squad',lambda:q.squad(q1,qnear,.37)),
                     ('squad4',lambda:quaternion.squad4(q,q1,qnear,q,.37)),
                     ('to_matrix',q.toMatrix),('from_matrix',lambda:quaternion.quat_from_matrix(qm))]:
        add('quaternion_'+label,fn)
    translation=matrix.Matrix(4).translate(V(2,-3,.5))
    rotation=matrix.Matrix(3).rotate(V(0,0,1),37)
    base4=matrix.Matrix(4);base3=matrix.Matrix(3);axis=V3.normalize()
    add('matrix4_translate',lambda:base4.translate(V3))
    add('matrix3_rotate',lambda:base3.rotate(axis,37))
    add('vector3_affine_transform',lambda:vector.transform(3,V3.vector,translation.matrix))
    controls=[V(0,0,0),V(.2,2,.5),V(.8,-2,-.5),V(1,0,0)]
    quadratic=controls[:3]
    add('bezier_quadratic',lambda:bezier.quadraticBezierPoint(.37,*quadratic))
    add('bezier_cubic',lambda:bezier.cubicBezierPoint(.37,*controls))
    add('bezier_cubic_scalar',lambda:bezier.cubicBezierPoint(.37,0,2,-2,1))
    for tolerance in (.1,.01,1e-4,1e-8,1e-15):
        p=path(tolerance)
        add(f'bezier_subdivide_{tolerance:g}',lambda p=p:p.findDrawingPoints(0),
            'adaptive depth sweep',tolerance_distance=tolerance,segments=1)
    p=path(.001,8)
    add('bezier_path_8_segments',p.getDrawingPoints,'multi-segment',segments=8)
    for l,m in [(0,0),(2,1),(8,3),(12,6)]:
        add(f'sh_basis_{l}_{m}',lambda l=l,m=m:sh.SPH(l,m,.73,1.27),degree=l,order=m)
    # Fixed midpoint lattice; no RNG and no environmental asset dependency.
    for count in (64,256,1024):
        for bands in (1,3,5):
            samples=[];colors=[]
            for i in range(count):
                z=1-2*(i+.5)/count;phi=(i*.6180339887498949%1)*2*math.pi
                theta=math.acos(z);d=V(math.sqrt(1-z*z)*math.cos(phi),math.sqrt(1-z*z)*math.sin(phi),z)
                s=sh.SPHSample(theta,phi,d,bands*bands)
                s.values[:]=[sh.SPH(l,m,theta,phi) for l in range(bands) for m in range(-l,l+1)]
                samples.append(s);colors.append([1+4*max(d.vector[0],0)**8,.5+max(z,0),.25])
            add(f'sh_project_n{count}_b{bands}',lambda samples=samples,colors=colors:sh.project_radiance(samples,colors),
                'projection scaling',samples=count,bands=bands,coefficients=bands*bands)
    rgb=[[.05*(i+1),.025*(-1)**i,.1/(i+1)] for i in range(9)]
    scalar=[row[0] for row in rgb]
    add('sh_convolve_l2',lambda:sh.convolve_diffuse(rgb))
    add('sh_rotate_l2_scalar',lambda:sh.rotate_coefficients(scalar,q))
    add('sh_rotate_l2_rgb',lambda:sh.rotate_coefficients(rgb,q))
    normal=V(0,0,1)
    add('sh_reconstruct_l2',lambda:sh.reconstruct(rgb,normal))
    p=plane.Plane();p.fromCoeffs(1,2,3,-4)
    add('plane_normalize',p.normalize)
    point=V(1,2,3,1)
    add('plane_dot',lambda:p.dot(point))
    construction=[V(0,0,2),V(1,0,2),V(0,1,2)]
    def plane_points():
        p=plane.Plane();p.fromPoints(*construction);return p
    add('plane_from_points',plane_points)
    r=ray.Ray(V(1,2,3),V(2,1,-1))
    add('ray_duplicate',r.duplicate)
    # Mutable transforms operate on fresh duplicates; timings include that cost.
    def ray_transform(method,arg):
        fresh=r.duplicate();getattr(fresh,method)(arg);return fresh
    add('ray_translate',lambda:ray_transform('translate',translation))
    add('ray_matrix_rotate',lambda:ray_transform('roateUsingMatrix',rotation))
    add('ray_quaternion_rotate',lambda:ray_transform('rotateUsingQuaternion',q))
    return calls,details


V3,V3B=V(1.25,-2.5,3.75),V(-.5,2,1)


def measure(fn,trials,target):
    fn() # warm, validates the executable path before calibration
    timer=timeit.Timer(fn,timer=time.perf_counter)
    count=1
    while True:
        elapsed=timer.timeit(count)
        if elapsed>=target:break
        count*=2
    samples=[timer.timeit(count)/count for _ in range(trials)]
    median=statistics.median(samples)
    return {'operations_per_trial':count,'trials':trials,'seconds_per_operation':samples,
            'median_us':median*1e6,'min_us':min(samples)*1e6,'max_us':max(samples)*1e6,
            'mad_us':statistics.median(abs(s-median) for s in samples)*1e6}


def profile(fn,count):
    profiler=cProfile.Profile();profiler.enable()
    for _ in range(count):fn()
    profiler.disable();stats=pstats.Stats(profiler)
    rows=[]
    for (file,line,name),(primitive,total,self_time,cumulative,callers) in stats.stats.items():
        rows.append({'function':f'{file}:{line}:{name}','primitive_calls':primitive,'calls':total,
                     'self_seconds':self_time,'cumulative_seconds':cumulative})
    rows.sort(key=lambda row:row['cumulative_seconds'],reverse=True)
    tracemalloc.start();fn();current,peak=tracemalloc.get_traced_memory();tracemalloc.stop()
    return {'operations':count,'top_cumulative':rows[:20],'single_call_peak_traced_bytes':peak,
            'single_call_retained_traced_bytes':current}


def historical(trials,target):
    text=subprocess.check_output(['git','show',HISTORICAL+':gem/matrix.py'],cwd=ROOT,text=True)
    ns={'__name__':'historical_matrix'};exec(compile(text,'historical_matrix.py','exec'),ns)
    comparison={}
    for n in (3,4):
        rows=[[float(2 if i==j else 1 if j==i+1 else 0) for j in range(n)] for i in range(n)]
        old=ns['inverse'+str(n)];new=getattr(matrix,'inverse'+str(n))
        oldobj=ns['Matrix'](n,rows);newobj=matrix.Matrix(n,rows)
        timings={name:measure(fn,trials,target) for name,fn in
                 [('historical_raw',lambda:old(rows)),('current_raw',lambda:new(rows)),
                  ('historical_wrapper',oldobj.inverse),('current_wrapper',newobj.inverse)]}
        checks=[]
        for scale in (1e-300,1e-200,1,1e200,1e300):
            scaled=[[v*scale for v in row] for row in rows]
            inv=new(scaled)
            expected=[[((-1)**(j-i)/2**(j-i+1))/scale if j>=i else 0 for j in range(n)] for i in range(n)]
            error=max(abs((inv[i][j]-expected[i][j])/expected[i][j]) for i in range(n) for j in range(i,n))
            assert error<1e-14
            orders=[matrix.matrix_multiply(scaled,inv),matrix.matrix_multiply(inv,scaled)]
            residual=max(abs(p[i][j]-(i==j)) for p in orders for i in range(n) for j in range(n))
            assert residual<1e-14
            try:
                previous=old(scaled)
                old_error=max(abs((previous[i][j]-expected[i][j])/expected[i][j]) for i in range(n) for j in range(i,n))
                old_status='finite accurate' if math.isfinite(old_error) and old_error<1e-14 else 'inaccurate/nonfinite'
            except (ZeroDivisionError,OverflowError):old_status='arithmetic exception'
            checks.append({'scale':scale,'current_max_relative_error':error,'current_two_order_identity_residual':residual,'historical_status':old_status})
        comparison[str(n)]={'timings':timings,'extreme_known_answers':checks,
            'raw_slowdown':timings['current_raw']['median_us']/timings['historical_raw']['median_us'],
            'wrapper_slowdown':timings['current_wrapper']['median_us']/timings['historical_wrapper']['median_us']}
    return {'commit':HISTORICAL,'results':comparison}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--trials',type=int,default=7)
    ap.add_argument('--target-seconds',type=float,default=.02)
    ap.add_argument('--historical',action='store_true',help='requires historical git object')
    ap.add_argument('--profiles',action='store_true')
    args=ap.parse_args()
    if args.trials<3 or not math.isfinite(args.target_seconds) or args.target_seconds<=0:ap.error('require >=3 trials and finite positive target')
    start=time.perf_counter();calls,details=prepare();setup=time.perf_counter()-start
    fingerprint=hashlib.sha256()
    for file in sorted((ROOT/'gem').rglob('*.py')):
        fingerprint.update(file.relative_to(ROOT).as_posix().encode());fingerprint.update(file.read_bytes())
    results={name:{**measure(fn,args.trials,args.target_seconds),**details[name]} for name,fn in calls.items()}
    for name,fn in calls.items():
        if name.startswith('bezier_subdivide'):
            results[name]['output_points']=len(fn())
    cold=[]
    code="import time; t=time.perf_counter(); from gem import vector,matrix,quaternion,bezier,spherical_harmonics,plane,ray; print(time.perf_counter()-t)"
    for _ in range(args.trials):
        cold.append(float(subprocess.check_output([sys.executable,'-c',code],cwd=ROOT,text=True)))
    cpu='unavailable'
    if Path('/proc/cpuinfo').exists():
        cpu=next((s.split(':',1)[1].strip() for s in Path('/proc/cpuinfo').read_text().splitlines() if s.startswith('model name')),cpu)
    data={'source':{'core_sha256':fingerprint.hexdigest(),'harness_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},'environment':{'python':sys.version,'implementation':platform.python_implementation(),
          'platform':platform.platform(),'architecture':platform.machine(),'cpu':cpu,
          'cpu_count':os.cpu_count(),'six':importlib.metadata.version('six')},
          'methodology':{'timer':'perf_counter wall time','gc':'timeit disables GC during trials',
          'trials':args.trials,'target_seconds':args.target_seconds,'setup_seconds':setup,
          'no_op_overhead':measure(lambda:None,args.trials,args.target_seconds),
          'cold_import_seconds':cold,'cold_import_median_seconds':statistics.median(cold),
          'cold_import_scope':'fresh interpreter internal imports; excludes process launch; warm OS cache',
          'allocation':'public results included; constructor cases separately measured; no subtraction',
          'seed':'no randomness; fixed source constants and deterministic lattice'},'results':results,
          'unsupported':{'ray_intersections':'no intersection method exists in supported Ray; not fabricated'}}
    if args.historical:data['historical_inverse']=historical(args.trials,args.target_seconds)
    if args.profiles:
        names=['matrix3_raw_inverse','matrix4_raw_inverse','matrix4_multiply','quaternion_squad4',
               'bezier_subdivide_1e-15','bezier_path_8_segments','sh_project_n1024_b5',
               'sh_rotate_l2_rgb','ray_quaternion_rotate']
        data['profiles']={name:profile(calls[name],max(1,min(2000,int(.1/(results[name]['median_us']/1e6))))) for name in names}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
    print(f'{len(results)} cases saved to {args.output}')


if __name__=='__main__':main()
