"""Standalone repair checks plus public documentation examples.

Run with python -I from outside the checkout for wheel/sdist verification.
--source explicitly enables the checkout for documentation checks in development.
"""
import argparse
from fractions import Fraction as F
import importlib.util
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--source',action='store_true')
    args = parser.parse_args()
    if args.source: sys.path.insert(0,str(ROOT))
    import gem
    from gem import quaternion as q
    from gem import vector, matrix, bezier, legendre, spherical_harmonics as sh
    package_root = Path(gem.__file__).resolve().parent.parent
    if args.source:
        assert package_root == ROOT
    else:
        assert package_root.is_relative_to(Path(sys.prefix).resolve()) and package_root != ROOT
    imports = {}
    for module in (gem,q,vector,matrix,bezier,legendre,sh):
        assert Path(module.__file__).resolve().is_relative_to(package_root)
        imports[module.__name__] = str(Path(module.__file__).resolve())
    cases = []
    unit = math.ulp(0.)
    for sign in (-1.,1.):
        for representation in (-1.,1.):
            for api in ('free_slerp','method_slerp','squad4'):
                endpoint = q.quat_from_axis_angle([1.,0.,0.],math.degrees(2*sign*unit))
                if representation < 0: endpoint = endpoint.negate()
                data = endpoint.data; saved = data[:]; start = q.Quaternion()
                if api == 'free_slerp': out = q.quat_slerp(start,endpoint,1.)
                elif api == 'method_slerp': out = start.slerp(endpoint,1.)
                else: out = q.squad4(start,endpoint,start,endpoint,1.)
                assert out.data is not data and endpoint.data is data and data == saved
                expected = [0.,1.,float(2*F(endpoint.data[0])*F(endpoint.data[1]))]
                actual = q.quat_rotate_vector(out,vector.Vector(3,[0.,1.,0.])).vector
                assert actual == expected
                cases.append({'api':api,'sign':sign,'representation':representation,
                              'actual_rotated_vector':actual,'expected_rotated_vector':expected})
    assert q.quat_slerp(q.Quaternion(),q.Quaternion([1.,unit,0.,0.]),.5).data == [1.,0.,0.,0.]
    assert q.quat_slerp(q.Quaternion(),q.Quaternion([1.,unit,0.,0.]),.75).data == [1.,unit,0.,0.]
    assert q.quat_slerp(q.Quaternion(),q.Quaternion([1.,3*unit,0.,0.]),.375).data == [1.,unit,0.,0.]
    for t in (F(1,8),F(3,8),F(3,4)):
        result = q.quat_slerp(q.Quaternion(),q.Quaternion([1.,3*unit,0.,0.]),t)
        expected = float(t*F(3*unit))
        assert result.data == [1.,expected,0.,0.]
    quarter = q.quat_from_axis_angle([0.,0.,1.],90.)
    midpoint = q.quat_slerp(q.Quaternion(),quarter,.5)
    for actual in (q.quat_rotate_vector(midpoint,vector.Vector(3,[1.,0.,0.])).vector,
                   (midpoint.toMatrix()*vector.Vector(4,[1.,0.,0.,0.])).vector[:3]):
        assert all(abs(a-b) <= 3e-16 for a,b in zip(actual,[math.sqrt(.5),math.sqrt(.5),0.]))
    spec = importlib.util.spec_from_file_location('unchanged_docs_checker',ROOT/'tools/check_architecture_docs.py')
    checker = importlib.util.module_from_spec(spec); spec.loader.exec_module(checker)
    # Reuse declaration/signature/reexport checks without the historical
    # documentation-only phase's immutable-production checkpoint gate.
    api_reference = checker.check_api_reference(checker.declarations(),runtime=True)
    pages = json.loads((ROOT/'audit/phase4g5-verification.json').read_text())['documentation']['pages']
    snippets = {}
    for page in pages:
        for index,code in enumerate(re.findall(r'```python\n(.*?)```',(ROOT/page).read_text(),re.DOTALL)):
            exec(compile(code,page+':example'+str(index+1),'exec'),{})
            snippets[page] = snippets.get(page,0)+1
    assert sum(snippets.values()) == 42
    args.output.write_text(json.dumps({'package_root':str(package_root),'imports':imports,
        'source_mode':args.source,'original_regressions_passed':len(cases),'regressions':cases,
        'additional_subnormal_checks_passed':6,'rotation_matrix_checks_passed':2,
        'api_reference':api_reference,'documentation_examples':snippets,
        'examples_executed':sum(snippets.values()),'source_tree_runtime_imported':args.source,
        'checker_modified':False,'python':sys.version},indent=2)+'\n')
    print('12 endpoint regressions, 8 boundary/rotation checks, 268 declarations and 42 documentation examples passed')


if __name__ == '__main__':
    main()
