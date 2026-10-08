"""Run from repository root: python benchmarks/baseline.py --output audit/benchmark-results.json."""
import argparse
import gc
import json
import platform
import statistics
import subprocess
import timeit
from pathlib import Path
from gem import matrix, quaternion, vector
from gem.experimental import bezier, legendre, sph


def run():
    a = vector.Vector(3,[1.0,2.0,3.0])
    b = vector.Vector(3,[4.0,5.0,6.0])
    m = matrix.Matrix(4,[[2,1,0,0],[0,3,1,0],[0,0,4,0],[2,3,4,1]])
    v = vector.Vector(4,[1,2,3,1])
    q = quaternion.quat_from_axis_angle(vector.Vector(3,[0,0,1]),45)
    cases = {
        'vector3_add': lambda: a+b,
        'vector3_dot': lambda: a.dot(b),
        'vector3_cross': lambda: vector.cross(a,b),
        'vector3_normalize': lambda: a.normalize(),
        'matrix4_raw_multiply': lambda: matrix.matrix_multiply(m.matrix,m.matrix),
        'matrix4_object_multiply': lambda: m*m,
        'matrix4_vector4_multiply': lambda: m*v,
        'matrix4_inverse': lambda: m.inverse(),
        'quaternion_multiply': lambda: q*q,
        'quaternion_rotate_vector3': lambda: quaternion.quat_rotate_vector(q,a),
        'quaternion_slerp': lambda: q.slerp(quaternion.Quaternion(),0.5),
        'quadratic_bezier_scalar': lambda: bezier.quadraticBezierPoint(0.5,0,1,2),
        'legendre_l2_m0': lambda: legendre.Legendre(2,0,0.5).run(),
        'spherical_harmonic_l2_m0': lambda: sph.SPH(2,0,0.5,0.5),
    }
    # Check benchmark inputs before timing; avoid timing broken paths as correct work.
    assert a.dot(b) == 32
    assert abs((m*m.inverse()).matrix[0][0]-1) < 1e-12
    assert abs(quaternion.quat_rotate_vector(q,a).magnitude()-a.magnitude()) < 1e-12
    results = {}
    for name, fn in cases.items():
        fn()  # warmup
        timer = timeit.Timer(fn)
        count, _ = timer.autorange()
        samples = [seconds/count*1e6 for seconds in timer.repeat(repeat=5,number=count)]
        results[name] = {'iterations_per_repeat':count,'repeats':5,
                         'samples_us':samples,'median_us':statistics.median(samples),
                         'min_us':min(samples),'max_us':max(samples)}
    return {'python':platform.python_version(),'implementation':platform.python_implementation(),
            'platform':platform.platform(),'processor':platform.processor(),
            'source_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
            'gc_during_timing':'disabled by timeit (restored afterward)',
            'units':'microseconds per operation','results':results}


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    result = run()
    payload = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(payload)
    print(payload)
