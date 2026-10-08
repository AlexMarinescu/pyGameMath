import math
import pytest
from gem import bezier, vector


def V(*values):
    return vector.Vector(len(values), list(values))


def path(points, squared_tolerance=0.01):
    result = bezier.BezierPath()
    result.setControlPoints(points)
    result.minimum_sqr_distance = squared_tolerance
    return result


def distance_to_segment(point, a, b):
    delta = [y-x for x,y in zip(a,b)]
    denominator = sum(x*x for x in delta)
    t = min(1,max(0,sum((p-x)*d for p,x,d in zip(point,a,delta))/denominator)) if denominator else 0
    return math.sqrt(sum((p-x-t*d)**2 for p,x,d in zip(point,a,delta)))


def reference_error(samples, reference):
    coords = [p.vector for p in samples]
    return max(min(distance_to_segment(reference(i/1000),a,b)
                   for a,b in zip(coords,coords[1:])) for i in range(1001))


@pytest.mark.parametrize('dimension', [2,3])
@pytest.mark.parametrize('quadratic', [False,True])
def test_independent_reference_order_error_and_tolerance(dimension, quadratic):
    # x=t; cubic y=3t(1-t)(1-2t), or elevated quadratic y=2t(1-t).
    heights = [0,2/3,2/3,0] if quadratic else [0,1,-1,0]
    controls = [V(i/3,y,*([2*i/3] if dimension==3 else [])) for i,y in enumerate(heights)]
    original = [p.vector[:] for p in controls]
    storage = [p.vector for p in controls]
    reference = lambda t: [t,2*t*(1-t) if quadratic else 3*t*(1-t)*(1-2*t)] + ([2*t] if dimension==3 else [])
    errors=[]; counts=[]
    for tolerance in [0.1,0.02,0.004]:
        samples = path(controls,tolerance*tolerance).findDrawingPoints(0)
        counts.append(len(samples))
        assert samples[0].vector == original[0] and samples[-1].vector == original[-1]
        assert all(a.vector[0]<b.vector[0] for a,b in zip(samples,samples[1:]))
        for sample in samples:
            assert sample.vector == pytest.approx(reference(sample.vector[0]),abs=1e-14)
            assert all(sample is not p and sample.vector is not p.vector for p in controls)
        error=reference_error(samples,reference)
        assert error <= tolerance + 1e-14
        errors.append(error)
    assert errors == sorted(errors,reverse=True)
    assert counts == sorted(counts)
    assert [p.vector for p in controls] == original
    assert all(p.vector is s for p,s in zip(controls,storage))


@pytest.mark.parametrize('dimension',[2,3])
def test_lines_coincidence_and_endpoint_loops(dimension):
    line=[V(i,*([0]*(dimension-1))) for i in range(4)]
    assert [p.vector for p in path(line).findDrawingPoints(0)] == [line[0].vector,line[-1].vector]
    point=V(*([2]*dimension))
    samples=path([point]*4).findDrawingPoints(0)
    assert len(samples)==2 and samples[0].vector==samples[1].vector==point.vector
    assert samples[0] is not samples[1]
    loop=path([V(*([0]*dimension)),V(1,*([1]*(dimension-1))),
               V(-1,*([1]*(dimension-1))),V(*([0]*dimension))]).findDrawingPoints(0)
    assert len(loop)>2
    assert loop[0].vector==loop[-1].vector==[0]*dimension


def test_collinear_overshoot_is_not_discarded():
    samples=path([V(0,0),V(4,0),V(4,0),V(1,0)],1e-6).findDrawingPoints(0)
    assert max(p.vector[0] for p in samples)>3


def test_connected_segments_nested_shape_and_repeated_calls():
    controls=[V(i,0) for i in range(7)]
    p=path(controls)
    assert [[q.vector for q in group] for group in p.getDrawingPoints()] == [[[0,0],[3,0]],[[6,0]]]
    assert [[q.vector for q in group] for group in p.getDrawingPoints()] == [[[0,0],[3,0]],[[6,0]]]
    assert p.getControlPoints() is controls
    assert path([]).getDrawingPoints()==[]


def test_scale_and_squared_tolerance_semantics():
    controls=[V(0,0),V(1/3,1),V(2/3,1),V(1,0)]
    small=path(controls,0.01).findDrawingPoints(0)
    scaled=[V(*(10*x for x in p.vector)) for p in controls]
    assert len(path(scaled,0.01).findDrawingPoints(0))>len(small)
    proportional=path(scaled,1).findDrawingPoints(0)
    assert len(proportional)==len(small)
    for a,b in zip(small,proportional):
        assert b.vector==pytest.approx([10*x for x in a.vector])


def test_subinterval_insertion_count_and_parameter_order():
    p=path([V(0,0),V(1/3,1),V(2/3,-1),V(1,0)],1e-5)
    left,right=V(.25,.28125),V(.75,-.28125)
    result=[left,right]
    added=p.findDrawingPointsAdded(0,.25,.75,result,1)
    assert added==len(result)-2 and added>0
    assert result[0] is left and result[-1] is right
    assert all(a.vector[0]<b.vector[0] for a,b in zip(result,result[1:]))
    for q in result:
        t=q.vector[0]
        assert q.vector[1]==pytest.approx(3*t*(1-t)*(1-2*t))


def test_actual_depth_limit_best_effort():
    p=path([V(0,0),V(1/3,1),V(2/3,1),V(1,0)],1e-300)
    samples=p.findDrawingPoints(0)
    assert len(samples)==65537
    assert samples[0].vector==[0,0] and samples[-1].vector==[1,0]
    assert all(a.vector[0]<b.vector[0] for a,b in zip(samples,samples[1:]))
    # At depth 16, the first chord's midpoint misses the parabola by much
    # more than the requested 1e-150 distance: capped output is best effort.
    a,b=samples[:2]; t=(a.vector[0]+b.vector[0])/2
    assert distance_to_segment([t,3*t*(1-t)],a.vector,b.vector)>1e-150


@pytest.mark.parametrize('points', [[V(0,0)]*3,[V(0,0)]*5,[V(0,0),V(1,2,3),V(2,0),V(3,0)]])
def test_malformed_sampling(points):
    with pytest.raises(ValueError):path(points).findDrawingPoints(0)


@pytest.mark.parametrize('tolerance',[0,-1,float('inf'),float('nan')])
def test_invalid_tolerances(tolerance):
    with pytest.raises(ValueError):path([0,1,2,3],tolerance).findDrawingPoints(0)


def test_builders_known_controls_ownership_append_and_rebuild():
    controls=[V(0,0),V(1,0),V(2,0)]
    saved=[p.vector[:] for p in controls]; storage=[p.vector for p in controls]
    p=bezier.BezierPath()
    assert p.interpolate(controls,.5) is None
    expected=[[0,0],[.5,0],[.5,0],[1,0],[1.5,0],[1.5,0],[2,0]]
    assert [q.vector for q in p.controlPoints]==expected
    assert p.curveCount==2
    assert p.interpolate(controls,.5) is None
    assert [q.vector for q in p.controlPoints]==expected+expected
    for _ in range(2):
        assert p.samplePoints(controls,.01,1,.5) is None
        assert [q.vector for q in p.controlPoints]==expected
        assert p.curveCount==2
    assert [p.vector for p in controls]==saved
    assert all(p.vector is s for p,s in zip(controls,storage))
    assert all(q is not v and q.vector is not v.vector for q in p.controlPoints for v in controls)
    old=p.controlPoints
    p.samplePoints([],0,1,.5);p.interpolate([controls[0]],.5)
    assert p.controlPoints is old


def test_thinning_and_coincident_builders():
    p=bezier.BezierPath()
    p.samplePoints([V(x,0) for x in [0,.1,.2,1]],.25,1,.5)
    assert p.curveCount==1
    assert [q.vector for q in p.controlPoints]==[[0,0],[.5,0],[.5,0],[1,0]]
    p.samplePoints([V(x,0) for x in [0,.1,.2,1]],.25,.5,.5)
    assert p.curveCount==2 # retain .2 before next gap exceeds max
    p.samplePoints([V(2,3)]*3,0,1,.5)
    assert all(q.vector==[2,3] for q in p.controlPoints)
    assert p.findDrawingPoints(0)[0].vector==[2,3]


def test_compatibility_class_is_core():
    from gem.experimental.bezier import BezierPath
    from gem.experimental._bezier_legacy import BezierPath as private
    assert BezierPath is private is bezier.BezierPath


@pytest.mark.parametrize('values,minimum,maximum,retained', [
    ([0,.5,1],.25,.25,[0,.5,1]),
    ([0,.01,.02,1],0,1,[0,.01,.02,1]),
    ([0,.01,.02,.03,1],.25,1,[0,1]),
    ([0,10,20],1,4,[0,10,20]),
    ([0,20],1,4,[0,20]),
])
def test_threshold_boundaries_sparse_dense_and_gap_limit(values,minimum,maximum,retained):
    p=bezier.BezierPath()
    p.samplePoints([V(x,0) for x in values],minimum,maximum,0)
    assert [q.vector[0] for q in p.controlPoints[::3]]==retained
    assert p.controlPoints[0].vector==[values[0],0]
    assert p.controlPoints[-1].vector==[values[-1],0]
    if values[-1]>=20:
        assert max(b-a for a,b in zip(retained,retained[1:]))**2>maximum


@pytest.mark.parametrize('minimum,maximum',[(-1,1),(2,1),(0,0),(0,float('inf')),(float('nan'),1)])
def test_invalid_thinning_preserves_existing_path(minimum,maximum):
    p=path([0,1,2,3]); old=p.controlPoints
    with pytest.raises(ValueError):p.samplePoints([V(0,0),V(1,0)],minimum,maximum,.5)
    assert p.controlPoints is old


def test_interpolate_does_not_extend_caller_control_list():
    controls=[V(x,0) for x in range(4)]
    p=path(controls)
    p.interpolate([V(4,0),V(5,0)],.5)
    assert len(controls)==4
    assert p.controlPoints[:4]==controls and len(p.controlPoints)==8
    assert all(a is b for a,b in zip(p.controlPoints,controls))
