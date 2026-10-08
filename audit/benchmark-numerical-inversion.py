"""Compare final inverses with Phase 2F-1; run from the repository root."""
import json
import math
import statistics
import subprocess
import sys
import time
import timeit
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from gem import matrix

BASE='2cbd899a47546b8c40f8244f57799969c86b0d87'
namespace={'__name__':'baseline_matrix'}
exec(compile(subprocess.check_output(['git','show',BASE+':gem/matrix.py'],text=True),'<baseline>', 'exec'),namespace)
result={}
for size in [3,4]:
    rows=[[float(2 if i==j else 1 if j==i+1 else 0) for j in range(size)] for i in range(size)]
    old=namespace['inverse'+str(size)];new=getattr(matrix,'inverse'+str(size))
    old_wrapper=namespace['Matrix'](size,[r[:] for r in rows])
    new_wrapper=matrix.Matrix(size,[r[:] for r in rows])
    calls={'baseline_raw':lambda:old(rows),'final_raw':lambda:new(rows),
           'baseline_wrapper':old_wrapper.inverse,'final_wrapper':new_wrapper.inverse}
    timings={label:statistics.median(timeit.repeat(fn,number=3000,repeat=5,timer=time.process_time))/3000*1e6 for label,fn in calls.items()}
    errors=[]
    for scale in [1e-300,1e-200,1,1e200,1e300]:
        inv=new([[c*scale for c in row] for row in rows])
        expected=[[((-1)**(j-i)/2**(j-i+1))/scale if j>=i else 0 for j in range(size)] for i in range(size)]
        error=max(abs((inv[i][j]-expected[i][j])/expected[i][j]) for i in range(size) for j in range(i,size))
        errors.append({'scale':scale,'max_relative_error':error})
    result[str(size)]={'process_cpu_microseconds_median_5x3000':timings,'known_answer_errors':errors}
print(json.dumps({'base_commit':BASE,'python':sys.version.split()[0],'results':result},indent=2))
