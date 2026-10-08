"""Measure SLERP accuracy/cost against the Phase 2E-4 base (run from repo root)."""
import json
import math
import statistics
import subprocess
import timeit
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gem import quaternion as q

BASE = '4eafc86731aac609a6bc81113207c0e5b9444308'
source = subprocess.check_output(['git', 'show', BASE + ':gem/quaternion.py'], text=True)
old_namespace = {'__name__': 'interpolation_baseline'}
exec(compile(source, '<baseline>', 'exec'), old_namespace)
old_slerp = old_namespace['quat_slerp']
old_type = old_namespace['Quaternion']

measurements = []
for degrees in [0.001, 0.1, 1, 3, math.degrees(2*math.acos(0.999))-1e-8]:
    alpha = math.radians(degrees)/2
    item = {'spatial_separation_degrees': degrees}
    for name, fn, cls in [('baseline', old_slerp, old_type), ('corrected', q.quat_slerp, q.Quaternion)]:
        a = cls(); b = cls([math.cos(alpha), 0, 0, math.sin(alpha)])
        norm_error = angle_error = component_error = 0.0
        for i in range(1001):
            t = i/1000
            out = fn(a, b, t).data
            norm_error = max(norm_error, abs(math.hypot(out[0], out[3])-1))
            angle_error = max(angle_error, 2*abs(math.atan2(out[3], out[0])-t*alpha))
            component_error = max(component_error, abs(out[0]-math.cos(t*alpha)), abs(out[3]-math.sin(t*alpha)))
        item[name] = {'max_norm_error': norm_error, 'max_spatial_angle_error_radians': angle_error,
                      'max_component_error': component_error}
    measurements.append(item)

cost = {}
for name, fn, cls in [('baseline', old_slerp, old_type), ('corrected', q.quat_slerp, q.Quaternion)]:
    a = cls(); b = cls([math.cos(math.pi/360), 0, 0, math.sin(math.pi/360)])
    cost[name] = statistics.median(timeit.repeat(lambda: fn(a, b, .37), number=50000, repeat=7))/50000*1e6
print(json.dumps({'base_commit': BASE, 'samples_per_separation': 1001,
                  'measurements': measurements, 'microseconds_per_call_median_7x50000': cost}, indent=2))
