"""Check preserved finite results and untouched public code against post-2G."""
import ast
import hashlib
import json
from pathlib import Path
import random
import struct
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from gem import matrix

BASE = '89cd4b97784d32625de65b8f465cfcbf6e102943'


def main():
    before = subprocess.check_output(['git', 'show', BASE+':gem/matrix.py'], cwd=ROOT, text=True)
    after = (ROOT/'gem/matrix.py').read_text()
    ns = {'__name__': 'baseline_matrix'};exec(compile(before,'baseline_matrix.py','exec'),ns)
    before_nodes = {n.name:n for n in ast.parse(before).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    after_nodes = {n.name:n for n in ast.parse(after).body if isinstance(n,(ast.FunctionDef,ast.ClassDef))}
    changed = {'_inverse3_cofactor','_inverse4_cofactor','_scaled_cofactor_inverse'}
    unchanged = [name for name in before_nodes if name not in changed]
    for name in unchanged:
        assert ast.dump(before_nodes[name]) == ast.dump(after_nodes[name]), name
    count = 0
    for size in (3,4):
        for seed in range(50):
            rng = random.Random(seed)
            original = [[3*(i==j)+rng.uniform(-.25,.25) for j in range(size)] for i in range(size)]
            for scale in (1e-300,1e-200,1.,1e200,1e300):
                rows = [[v*scale for v in row] for row in original]
                old = ns['inverse'+str(size)](rows)
                new = getattr(matrix,'inverse'+str(size))(rows)
                assert all(struct.pack('!d',a)==struct.pack('!d',b)
                           for ar,br in zip(new,old) for a,b in zip(ar,br))
                count += 1
    result = {'baseline_commit':BASE,'bit_identical_finite_inverses':count,
              'unchanged_ast_definitions':unchanged,'current_source_sha256':hashlib.sha256(after.encode()).hexdigest(),
              'scope':'differential preservation, separate from independent Fraction/integer tests'}
    Path(sys.argv[1]).write_text(json.dumps(result,indent=2)+'\n')
    print(f'{count} finite inverses bit-identical; {len(unchanged)} public/unrelated definitions unchanged')


if __name__ == '__main__':main()
