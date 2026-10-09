"""Independent geometry audit: exact areas, affine coordinates and Rodrigues.

No gem rotation, cross product or normalization supplies an expected answer.
Local line/plane calculations exercise existing operations, not a new core API.
"""
from decimal import Decimal, localcontext
from fractions import Fraction
import math
import random

import pytest

from gem import matrix, plane, quaternion, ray, vector


def V(values):
    return vector.Vector(len(values), list(values))


def F(x):
    return Fraction(x)


def dot(a, b):
    return sum((F(x)*F(y) for x,y in zip(a,b)), Fraction())


def subtract(a,b):
    return [F(x)-F(y) for x,y in zip(a,b)]


def cross(a,b):
    return [a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]]


def decimal_unit(values):
    with localcontext() as ctx:
        ctx.prec = 110
        d=[Decimal(F(x).numerator)/Decimal(F(x).denominator) for x in values]
        length=sum(x*x for x in d).sqrt()
        return [float(x/length) for x in d],float(length)


def polygon_area_vector(points):
    # Exact translated triangle fan, independent of Newell's sum/difference form.
    anchor=points[0]
    area=[Fraction()]*3
    for i in range(1,len(points)-1):
        triangle=cross(subtract(points[i],anchor),subtract(points[i+1],anchor))
        area=[x+y for x,y in zip(area,triangle)]
    return area


def signed_area(a,b,c):
    return (F(b[0])-F(a[0]))*(F(c[1])-F(a[1]))-(F(b[1])-F(a[1]))*(F(c[0])-F(a[0]))


def barycentric_reference(p,a,b,c):
    # Independent 2D signed-area ratios for the affine embeddings tested below.
    denominator=signed_area(a,b,c)
    return [signed_area(p,b,c)/denominator,signed_area(a,p,c)/denominator,
            signed_area(a,b,p)/denominator]


def coefficients(p):
    return [p.a,p.b,p.c,p.d]


def assert_preserved(objects,storages,values):
    for obj,storage,value in zip(objects,storages,values):
        assert obj.vector is storage and obj.vector == value


@pytest.mark.parametrize('seed',range(20))
def test_plane_exact_graph_incidence_homogeneous_distance_and_ownership(seed):
    rng=random.Random(440001+seed)
    A,B,C=[rng.randint(-6,6) for _ in range(3)]
    points=[[0,0,C],[2,0,2*A+C],[0,3,3*B+C]]
    inputs=list(map(V,points));storages=[v.vector for v in inputs]
    expected_normal,_=decimal_unit([-A,-B,1])
    p=plane.Plane();assert p.fromPoints(*inputs) is None
    expected=expected_normal+[-C*expected_normal[2]]
    assert coefficients(p)==pytest.approx(expected,rel=3e-15,abs=3e-15)
    assert p.normal.vector==pytest.approx(expected_normal,rel=3e-15,abs=3e-15)
    for _ in range(4):
        x,y=rng.randint(-10,10),rng.randint(-10,10)
        xyz=[x,y,A*x+B*y+C]
        assert dot([-A,-B,1,-C],xyz+[1])==0
        assert p.dot(V(xyz+[1]))==pytest.approx(0,abs=1e-13)
        for offset in [-2,0,3]:
            off=xyz[:];off[2]+=offset
            assert p.dot(V(off+[1]))==pytest.approx(offset*expected_normal[2],abs=1e-13)
    assert_preserved(inputs,storages,points)
    assert all(p.normal.vector is not storage for storage in storages)


@pytest.mark.parametrize('scale',[-8.,-.125,.125,1.,8.])
@pytest.mark.parametrize('w',[-2.,0.,1.,3.])
def test_plane_scaled_homogeneous_equation_clone_flip_and_normalize(scale,w):
    raw=[2*scale,3*scale,6*scale,-26*scale]
    p=plane.Plane();assert p.fromCoeffs(*raw) is None
    normal=p.normal;storage=normal.vector
    point=[1,2,3,w]
    assert p.dot(V(point))==float(dot(raw,point))
    q=p.normalize();expected=[x/(7*abs(scale)) for x in raw]
    assert coefficients(q)==pytest.approx(expected,rel=2e-15,abs=2e-15)
    assert q.normal.vector==pytest.approx(expected[:3],rel=2e-15,abs=2e-15)
    assert coefficients(p)==raw and p.normal is normal and normal.vector is storage
    assert q.normal is not normal and q.normal.vector is not storage
    for result,sign in [(p.clone(),1),(p.flip(),-1)]:
        assert coefficients(result)==[sign*x for x in raw]
        assert result.normal.vector==[sign*x for x in raw[:3]]
        assert result.normal is not normal and result.normal.vector is not storage
        result.normal.vector[0]=999
        assert p.normal.vector==raw[:3]
    receiver=plane.Plane();receiver.fromCoeffs(*raw)
    assert receiver.i_normalize() is receiver
    assert coefficients(receiver)==coefficients(q)
    assert receiver.i_flip() is receiver
    assert coefficients(receiver)==[-x for x in coefficients(q)]


@pytest.mark.parametrize('offset',[-3,0,7])
@pytest.mark.parametrize('scale',[-2,1,3])
def test_point_location_uses_explicit_plane_and_exact_sign(offset,scale):
    p=plane.Plane();p.fromCoeffs(2*scale,3*scale,6*scale,-26*scale)
    receiver=plane.Plane();receiver.fromCoeffs(1,0,0,999)
    point=[1,2,3+offset]
    exact=dot([2*scale,3*scale,6*scale,-26*scale],point+[1])
    assert receiver.point_location(p,tuple(point))==(1 if exact>0 else -1 if exact<0 else 0)
    assert coefficients(receiver)==[1,0,0,999]


@pytest.mark.parametrize('shape',[
    [[0,0,0],[2,0,0],[2,3,0],[0,3,0]],
    [[0,0,0],[2,0,2],[2,3,5],[0,3,3]],
    [[0,0,0],[3,0,0],[3,1,1],[0,1,2]], # nonplanar: area vector, not least squares
])
@pytest.mark.parametrize('reverse',[False,True])
@pytest.mark.parametrize('closed',[False,True])
def test_polygon_exact_area_winding_wrapping_and_weighted_offset(shape,reverse,closed):
    points=[p[:] for p in (list(reversed(shape)) if reverse else shape)]
    expected,_=decimal_unit(polygon_area_vector(points))
    if closed:points.append(points[0][:])
    vs=list(map(V,points));storages=[v.vector for v in vs]
    helper=plane.Plane();normal=helper.bestFitNormal(vs)
    assert normal.vector==pytest.approx(expected,rel=3e-15,abs=3e-15)
    expected_D=float(sum((dot(p,normal.vector) for p in points),Fraction())/len(points))
    assert helper.bestFitD(vs,normal)==pytest.approx(expected_D,rel=3e-15,abs=3e-15)
    for shift in range(len(shape)):
        shifted=points[:len(shape)][shift:]+points[:len(shape)][:shift]
        assert helper.bestFitNormal(list(map(V,shifted))).vector==pytest.approx(expected,abs=3e-15)
    assert_preserved(vs,storages,points)
    assert coefficients(helper)==[0,0,0,0]
    normal.vector[0]=999;assert_preserved(vs,storages,points)


@pytest.mark.parametrize('seed',range(20))
@pytest.mark.parametrize('dimension',[2,3,4])
def test_exact_affine_barycentric_reconstruction_orientation_and_ownership(seed,dimension):
    rng=random.Random(440100+seed)
    while True:
        points=[[rng.randint(-8,8),rng.randint(-8,8)] for _ in range(3)]
        if abs(signed_area(*points))>=8:break
    def embed(p):
        return p+[2*p[0]-p[1]+3,-p[0]+3*p[1]-2][:dimension-2]
    weights=[Fraction(rng.randint(-4,8),4),Fraction(rng.randint(-4,8),4)]
    weights.append(1-sum(weights))
    point=[sum(w*F(v[i]) for w,v in zip(weights,points)) for i in range(2)]
    expected=barycentric_reference(point,*points)
    assert expected==weights
    values=[list(map(float,embed(p))) for p in points+[point]]
    vs=list(map(V,values));storages=[v.vector for v in vs]
    a,b,c,p=vs
    result=p.barycentric(a,b,c)
    assert isinstance(result,list) and result==pytest.approx(list(map(float,expected)),rel=3e-12,abs=3e-12)
    assert sum(result)==pytest.approx(1,abs=3e-15)
    rebuilt=[sum(weight*vertex[i] for weight,vertex in zip(result,values[:3])) for i in range(dimension)]
    assert rebuilt==pytest.approx(list(map(float,values[3])),rel=5e-12,abs=5e-12)
    assert p.barycentric(a,c,b)==pytest.approx([float(weights[0]),float(weights[2]),float(weights[1])],abs=3e-12)
    assert p.barycentric(a,b,c) is not result
    result[0]=999;assert_preserved(vs,storages,values)


@pytest.mark.parametrize('seed',range(15))
def test_barycentric_off_plane_is_orthogonal_projection(seed):
    rng=random.Random(440200+seed)
    z=rng.randint(-4,4);h=rng.randint(-8,8)
    a,b,c=V([0,0,z]),V([3,0,z]),V([1,2,z])
    x,y=rng.randint(-4,4),rng.randint(-4,4)
    # XY signed areas independently give coordinates of orthogonal projection.
    expected=list(map(float,barycentric_reference([x,y],a.vector,b.vector,c.vector)))
    assert V([x,y,z+h]).barycentric(a,b,c)==pytest.approx(expected,abs=2e-15)


@pytest.mark.parametrize('normal',[[1.,0.,0.],[0.,-1.,0.],[2/7,3/7,6/7]])
@pytest.mark.parametrize('incident',[[2.,-3.,6.],[-.25,1.5,2.],[0.,0.,0.]])
def test_reflection_fraction_answer_involution_and_input_storage(normal,incident):
    i,n=V(incident),V(normal);storages=[i.vector,n.vector]
    product=dot(incident,normal)
    expected=[float(F(x)-2*product*F(y)) for x,y in zip(incident,normal)]
    out=vector.reflect(i,n)
    assert isinstance(out,vector.Vector) and out.vector==pytest.approx(expected,rel=3e-15,abs=3e-15)
    assert vector.reflect(out,n).vector==pytest.approx(incident,rel=3e-15,abs=3e-15)
    assert math.hypot(*out.vector)==pytest.approx(math.hypot(*incident),rel=3e-15,abs=3e-15)
    assert out is not i and out.vector is not i.vector and out.vector is not n.vector
    out.vector[0]=999;assert_preserved([i,n],storages,[incident,normal])


@pytest.mark.parametrize('seed',range(20))
def test_refraction_independent_snell_frame_tir_and_ownership(seed):
    rng=random.Random(440300+seed)
    phi=rng.uniform(-math.pi,math.pi);s,c=math.sin(phi),math.cos(phi)
    normal=[c,s,0.];tangent=[-s,c,0.]
    eta=rng.uniform(.5,2.);theta=rng.uniform(.1,.6)*math.asin(min(1,1/eta))
    incident=[math.sin(theta)*t-math.cos(theta)*n for t,n in zip(tangent,normal)]
    outgoing=math.asin(eta*math.sin(theta))
    expected=[math.sin(outgoing)*t-math.cos(outgoing)*n for t,n in zip(tangent,normal)]
    i,n=V(incident),V(normal);storages=[i.vector,n.vector]
    out=vector.refract(eta,i,n)
    assert out.vector==pytest.approx(expected,rel=3e-14,abs=3e-14)
    assert math.hypot(*out.vector)==pytest.approx(1,abs=3e-14)
    assert_preserved([i,n],storages,[incident,normal])
    # Independent transmitted tangential sine exceeds 1: no real refraction.
    beyond=[.9*t-math.sqrt(1-.9**2)*n for t,n in zip(tangent,normal)]
    tir=vector.refract(2.,V(beyond),n)
    assert tir.vector==[0.,0.,0.] and tir is not n


@pytest.mark.parametrize('dimension',[2,3,4])
@pytest.mark.parametrize('normal_component',[1e-9,1e-100])
@pytest.mark.defect('4G4-A01: equal-index grazing refraction loses nonzero normal component')
def test_equal_media_grazing_refraction_identity(dimension,normal_component):
    # Unit inputs to binary64 precision, as with ordinary sin/cos directions.
    incident=[math.sqrt(1-normal_component**2),-normal_component]+[0.]*(dimension-2)
    normal=[0.,1.]+[0.]*(dimension-2)
    assert math.hypot(*incident)==1.
    i,n=V(incident),V(normal);storages=[i.vector,n.vector]
    out=vector.refract(1.,i,n)
    assert_preserved([i,n],storages,[incident,normal])
    assert out is not i and out.vector is not i.vector
    # Snell with n1=n2: incidence equals transmission, including the tiny component.
    assert out.vector==pytest.approx(incident,rel=3e-15,abs=0)


@pytest.mark.parametrize('oblique',[False,True])
@pytest.mark.parametrize('reverse',[False,True])
@pytest.mark.parametrize('closed',[False,True])
@pytest.mark.defect('4G4-A02: translated Newell normal cancels for a nondegenerate polygon')
def test_translated_polygon_retains_exact_area_normal(oblique,reverse,closed):
    local=[[0.,0.,0.],[1.,0.,float(oblique)],[0.,1.,float(oblique)]]
    origin=float(2**52)
    points=[[origin+x for x in p] for p in local]
    assert all(subtract(p,[origin]*3)==list(map(F,q)) for p,q in zip(points,local))
    if reverse:points.reverse()
    if closed:points.append(points[0][:])
    area=polygon_area_vector(points)
    assert area==[F((-1 if oblique else 0)*(-1 if reverse else 1))]*2+[F(-1 if reverse else 1)]
    expected,_=decimal_unit(area)
    vs=list(map(V,points));storages=[v.vector for v in vs]
    # No overflow/underflow or nearly collinear triangle; exact area is nonzero.
    out=plane.Plane().bestFitNormal(vs)
    assert out.vector==pytest.approx(expected,rel=3e-15,abs=3e-15)
    assert_preserved(vs,storages,points)


def rodrigues(values,axis,angle):
    unit,_=decimal_unit(axis);s,c=math.sin(angle),math.cos(angle)
    product=sum(a*b for a,b in zip(unit,values));t=cross(unit,values)
    return [v*c+w*s+a*product*(1-c) for v,w,a in zip(values,t,unit)]


@pytest.mark.parametrize('axis',[[1,0,0],[0,1,0],[0,0,1],[2,-3,6]])
@pytest.mark.parametrize('angle',[0.,math.pi/2,-math.pi/2,math.pi,.73])
def test_ray_rodrigues_rigid_geometry_parameter_and_input_ownership(axis,angle):
    start=[1.,-2.,3.];raw=[2.,3.,6.];expected_dir,length=decimal_unit(raw)
    start_v,dir_v=V(start),V(raw);start_store,old_dir_store=start_v.vector,dir_v.vector
    r=ray.Ray(start_v,dir_v)
    assert r.start is start_v and r.dir is dir_v
    assert dir_v.vector is not old_dir_store and old_dir_store==raw
    assert r.distance==pytest.approx(length,rel=3e-15)
    assert r.dir.vector==pytest.approx(expected_dir,abs=3e-15)
    unit,_=decimal_unit(axis);half=angle/2
    qdata=[math.cos(half)]+[x*math.sin(half) for x in unit]
    q=quaternion.Quaternion(qdata[:]);qstore=q.data
    state=r.end;state.vector=[7.,8.,9.];endstore=state.vector
    assert r.rotateUsingQuaternion(q) is None
    expected_start=rodrigues(start,axis,angle);rotated_dir=rodrigues(expected_dir,axis,angle)
    assert r.start.vector==pytest.approx(expected_start,abs=1e-14)
    assert r.dir.vector==pytest.approx(rotated_dir,abs=1e-14)
    for parameter in [-3.,-0.,0.,.25,7.]:
        actual=r.start+r.dir*parameter
        expected=[p+d*parameter for p,d in zip(expected_start,rotated_dir)]
        assert actual.vector==pytest.approx(expected,abs=1e-13)
    assert r.distance==pytest.approx(length) and r.end is state and state.vector is endstore
    assert state.vector==[7.,8.,9.]
    assert q.data is qstore and q.data==qdata
    assert start_v.vector is start_store and start_v.vector==start
    assert dir_v.vector==pytest.approx(expected_dir,abs=3e-15)
    clone=r.duplicate()
    for original,duplicate in zip([r.start,r.dir,r.end],[clone.start,clone.dir,clone.end]):
        assert original is not duplicate and original.vector is not duplicate.vector
        assert original.vector==duplicate.vector
    assert clone.distance==r.distance


@pytest.mark.parametrize('scale',[1e-300,1.,1e300])
def test_ray_stable_constructor_direction_and_geometric_parameter(scale):
    raw=[2*scale,3*scale,6*scale];expected,length=decimal_unit(raw)
    start=V([1.,2.,3.]);direction=V(raw);r=ray.Ray(start,direction)
    assert r.distance==pytest.approx(length,rel=3e-15,abs=0)
    assert r.dir.vector==pytest.approx(expected,rel=3e-15,abs=0)
    assert (r.start+r.dir*7).vector==pytest.approx([3.,5.,9.],rel=3e-15)
    assert r.end.vector==[0.,0.,0.]


@pytest.mark.parametrize('order',['rotate-translate','translate-rotate'])
def test_literal_rigid_matrices_noncommute_and_preserve_incidence(order):
    rotation=[[0,1,0],[-1,0,0],[0,0,1]]
    rows=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[4,-5,2,1]]
    R,T=matrix.Matrix(3,[r[:] for r in rotation]),matrix.Matrix(4,[r[:] for r in rows])
    r=ray.Ray(V([1,2,3]),V([2,3,6]));dist=r.distance;end=r.end
    expected_start=[2,-4,5] if order=='rotate-translate' else [3,5,5]
    if order=='rotate-translate':r.roateUsingMatrix(R);r.translate(T)
    else:r.translate(T);r.roateUsingMatrix(R)
    assert r.start.vector==pytest.approx(expected_start,abs=2e-15)
    assert r.dir.vector==pytest.approx([-3/7,2/7,6/7],abs=2e-15)
    assert r.distance==dist and r.end is end and end.vector==[0.,0.,0.]
    # Rotate n=(2,3,6) to (-3,2,6), translate plane offset independently.
    moved_plane=plane.Plane();normal=[-3,2,6]
    translation=[4,-5,2] if order=='rotate-translate' else [5,4,2]
    d=-26-float(dot(normal,translation));moved_plane.fromCoeffs(*normal,d)
    assert moved_plane.dot(V(r.start.vector+[1]))==pytest.approx(0,abs=1e-13)
    for t in [-2,0,3]:
        p=r.start+r.dir*t
        assert moved_plane.dot(V(p.vector+[1]))==pytest.approx(7*t,abs=1e-13)
    assert R.matrix==rotation and T.matrix==rows


@pytest.mark.parametrize('slope',[0.,2**-30,-2**-30,1.])
@pytest.mark.parametrize('height',[0.,-3.,3.])
def test_tutorial_local_line_plane_relation_and_unbounded_parameters(slope,height):
    # There is no Ray intersection API: solve the tutorial equation locally.
    r=ray.Ray(V([-3,1,2+height]),V([1,0,slope]))
    p=plane.Plane();p.fromCoeffs(0,0,1,-2)
    value=p.dot(V(r.start.vector+[1]));denominator=p.normal.dot(r.dir)
    if slope==0:
        assert denominator==0 and value==height
        # Analytic coplanar/parallel distinction; do not change end state.
        assert (value==0)==(height==0)
    else:
        raw_parameter=-F(height)/F(slope)
        unit,length=decimal_unit([1,0,slope])
        parameter=-value/denominator
        expected_parameter=float(raw_parameter)*length
        assert parameter==pytest.approx(expected_parameter,rel=3e-15,abs=0)
        hit=(r.start+r.dir*parameter).vector
        assert hit==pytest.approx([-3+float(raw_parameter),1,2],rel=3e-15,abs=3e-15)
        assert p.dot(V(hit+[1]))==pytest.approx(0,abs=3e-15)
    assert r.end.vector==[0.,0.,0.] and r.distance==pytest.approx(math.hypot(1,slope))


@pytest.mark.parametrize('case',['ray','coeff-normal','three-points','polygon','mean','lookat-eye','lookat-up'])
def test_exact_degeneracy_retains_historical_errors(case):
    zero=V([0.,-0.,0.]);x=V([1,0,0]);p=plane.Plane()
    actions={'ray':lambda:ray.Ray(x,zero),'coeff-normal':lambda:p.normalize(),
             'three-points':lambda:p.fromPoints(zero,x,V([2,0,0])),
             'polygon':lambda:p.bestFitNormal([zero,x,zero]),
             'mean':lambda:p.bestFitD([],x),
             'lookat-eye':lambda:matrix.lookAt(x,x,V([0,1,0])),
             'lookat-up':lambda:matrix.lookAt(zero,x,x)}
    with pytest.raises(ZeroDivisionError):actions[case]()
    assert p.normal.vector==[0.,0.,0.] and coefficients(p)==[0,0,0,0]


def test_constructor_shared_vector_alias_is_established_mutation_not_copy():
    original=[0.,0.,5.];value=V(original);old=value.vector
    r=ray.Ray(value,value)
    assert r.start is r.dir is value
    assert r.start.vector==[0.,0.,1.] and value.vector is not old and old==original
    assert r.distance==5. and r.end.vector==[0.,0.,0.]
    duplicate=r.duplicate()
    assert duplicate.start is not duplicate.dir
    assert duplicate.start.vector is not duplicate.dir.vector
    assert duplicate.start.vector==duplicate.dir.vector==[0.,0.,1.]


@pytest.mark.parametrize('value',[[-3.,4.],[0.,-0.],[.25,-.5]])
def test_raw_2d_geometric_helpers_and_perpendicular_orientation(value):
    original=value[:]
    left,right=vector.lperp(value),vector.rperp(value)
    assert left.vector==[-value[1],value[0]]
    assert right.vector==[value[1],-value[0]]
    assert vector.toAngle(value)==math.atan2(value[1],value[0])
    assert dot(value,left.vector)==dot(value,right.vector)==0
    assert value==original and left.vector is not value and right.vector is not value


@pytest.mark.parametrize('scale',[-1e300,-1e-300,1e-300,1.,1e300])
def test_plane_normalization_finite_norm_scale_reference(scale):
    p=plane.Plane();raw=[2*scale,3*scale,6*scale,-26*scale];p.fromCoeffs(*raw)
    n=p.normalize();expected=[math.copysign(1,scale)*x/7 for x in [2,3,6,-26]]
    assert coefficients(n)==pytest.approx(expected,rel=3e-15,abs=0)
    assert coefficients(p)==raw and p.normal.vector==raw[:3]
    assert n.dot(V([1,2,3,1]))==pytest.approx(0,abs=2e-15)


@pytest.mark.parametrize('scale',[1e-150,1.,1e150])
def test_three_points_scaled_cross_with_representable_intermediates(scale):
    points=[[0,0,0],[scale,0,scale],[0,scale,scale]]
    expected,_=decimal_unit([-1,-1,1])
    p=plane.Plane();p.fromPoints(*map(V,points))
    assert p.normal.vector==pytest.approx(expected,rel=3e-15,abs=0)
    assert p.d==0
    for point in points:
        assert p.dot(V(point+[1]))==pytest.approx(0,abs=1e-15*scale)


@pytest.mark.parametrize('dimension',[2,4])
def test_axis_reflection_matching_dimensions_and_ownership(dimension):
    values=[2.,-3.]+[4.,5.][:dimension-2];normal=[0.,1.]+[0.]*(dimension-2)
    i,n=V(values),V(normal);out=vector.reflect(i,n)
    assert out.size==dimension and out.vector==[2.,3.]+values[2:]
    assert out.vector is not i.vector and out.vector is not n.vector
    assert i.vector==values and n.vector==normal


def test_lookat_nonaxis_known_frame_and_input_preservation():
    values=[[3.,4.,0.],[0.,0.,0.],[0.,0.,1.]]
    inputs=list(map(V,values));storages=[v.vector for v in inputs]
    view=matrix.lookAt(*inputs)
    expected=[[-.8,0,.6,0],[.6,0,.8,0],[0,1,0,0],[0,0,-5,1]]
    for row,answer in zip(view.matrix,expected):assert row==pytest.approx(answer,abs=2e-15)
    assert (view*V(values[0]+[1])).vector==pytest.approx([0,0,0,1],abs=2e-15)
    assert (view*V(values[1]+[1])).vector==pytest.approx([0,0,-5,1],abs=2e-15)
    assert_preserved(inputs,storages,values)


@pytest.mark.parametrize('transverse',[2**-30,1e-100,1e-300])
def test_structured_nearly_parallel_lookat_retains_nonzero_cross(transverse):
    view=matrix.lookAt(V([0,0,0]),V([1,0,0]),V([1,transverse,0]))
    expected=[[0,0,-1,0],[0,1,0,0],[1,0,0,0],[0,0,0,1]]
    for row,answer in zip(view.matrix,expected):assert row==pytest.approx(answer,abs=2e-15)
