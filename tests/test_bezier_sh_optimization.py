"""Independent polynomial, geometric and canonical-basis performance guards."""
import math
from fractions import Fraction
import pytest
from gem import bezier as b, spherical_harmonics as sh
from gem.vector import Vector
from gem.quaternion import Quaternion


def casteljau(values,t):
    values=list(map(Fraction,values));t=Fraction(t)
    while len(values)>1:
        values=[(1-t)*a+t*c for a,c in zip(values,values[1:])]
    return float(values[0])


@pytest.mark.parametrize('degree',[2,3])
@pytest.mark.parametrize('size',[0,1,2,3,4,8])
@pytest.mark.parametrize('t',[-.5,0.,.125,.5,1.,1.5])
def test_native_polynomial_exact_reference(degree,size,t):
    data=[[float((i+1)*(j+2)*(-1)**i) for j in range(size)] for i in range(degree+1)]
    controls=[Vector(size,p) for p in data]
    fn=b.quadraticBezierPoint if degree==2 else b.cubicBezierPoint
    result=fn(t,*controls)
    expected=[casteljau([p[j] for p in data],t) for j in range(size)]
    assert result.vector==expected
    assert result.size==size and all(result.vector is not p for p in data)
    assert all(obj.vector is p for obj,p in zip(controls,data))
    assert [obj.vector for obj in controls]==data


def test_polynomial_subclass_dispatch_and_mixed_dimensions():
    events=[]
    class Custom(Vector):
        def __mul__(self,value):
            events.append(value)
            return super().__mul__(value)
    controls=[Custom(2,[0.,0.]),Custom(2,[2.,4.]),Custom(2,[4.,0.])]
    assert b.quadraticBezierPoint(.5,*controls).vector==[2.,2.]
    assert events==[.25,.5,.25]
    # Preserve historical first-operand dimensionality when later controls are larger.
    assert b.quadraticBezierPoint(.5,Vector(2,[0.,0.]),Vector(3,[2.,4.,9.]),Vector(3,[4.,0.,8.])).vector==[2.,2.]


@pytest.mark.parametrize('polygon,expected',[
    ([(0.,0.,0.),(1.,3.,4.),(3.,-3.,-4.),(4.,0.,0.)],5.),
    ([(0.,0.),(-2.,0.),(6.,0.),(4.,0.)],2.),
    ([(0.,0.),(3.,4.),(-3.,-4.),(0.,0.)],5.),
    ([(0.,),(1.,),(2.,),(3.,)],0.),
])
def test_geometric_flatness_known_answers(polygon,expected):
    snapshot=polygon[:]
    assert b._flatness(polygon)==expected
    assert polygon==snapshot


@pytest.mark.parametrize('tolerance',[.03,.01,.003])
def test_adaptive_reference_error_and_parameter_order(tolerance):
    coords=[[0.,0.],[1/3,2.],[2/3,-2.],[1.,0.]]
    controls=[Vector(2,p) for p in coords]
    path=b.BezierPath();path.setControlPoints(controls);path.minimum_sqr_distance=tolerance*tolerance
    points=path.findDrawingPoints(0)
    assert points[0].vector==coords[0] and points[-1].vector==coords[-1]
    assert all(a.vector[0]<c.vector[0] for a,c in zip(points,points[1:]))
    for point in points:
        t=point.vector[0]
        assert point.vector[1]==pytest.approx(6*t*(1-t)*(1-2*t),abs=2e-14)
    def distance(p,a,c):
        dx,dy=c[0]-a[0],c[1]-a[1]
        u=max(0,min(1,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/(dx*dx+dy*dy)))
        return math.hypot(p[0]-a[0]-u*dx,p[1]-a[1]-u*dy)
    for i in range(201):
        t=i/200;p=(t,6*t*(1-t)*(1-2*t))
        assert min(distance(p,a.vector,c.vector) for a,c in zip(points,points[1:]))<=tolerance
    assert all(obj.vector is p for obj,p in zip(controls,coords))
    assert all(sample.vector is not p for sample in points for p in coords)


def test_depth_cap_and_exact_midpoint_reference():
    controls=[Vector(2,p) for p in [[0.,0.],[1/3,1.],[2/3,1.],[1.,0.]]]
    path=b.BezierPath();path.setControlPoints(controls);path.minimum_sqr_distance=1e-30
    points=path.findDrawingPoints(0)
    assert len(points)==2**16+1
    for index in (0,1,16384,32768,65535,65536):
        t=index/65536
        assert points[index].vector==pytest.approx([t,3*t*(1-t)],abs=3e-15)
    assert all(a.vector[0]<c.vector[0] for a,c in zip(points,points[1:]))


def test_degenerate_curve_and_shared_joins():
    controls=[Vector(3,[1.,2.,3.]) for _ in range(7)]
    path=b.BezierPath();path.setControlPoints(controls)
    points=path.getDrawingPoints()
    assert [len(segment) for segment in points]==[2,1]
    assert all(p.vector==[1,2,3] for segment in points for p in segment)
    assert points[0][0].vector is not points[0][1].vector
    assert all(p.vector is not c.vector for segment in points for p in segment for c in controls)


def basis_xyz(d):
    x,y,z=d;a=math.sqrt(3/(4*math.pi));c=math.sqrt(15/(4*math.pi))
    return [1/math.sqrt(4*math.pi),-a*y,a*z,-a*x,c*x*y,-c*y*z,
            math.sqrt(5/(16*math.pi))*(3*z*z-1),-c*x*z,math.sqrt(15/(16*math.pi))*(x*x-y*y)]


@pytest.mark.parametrize('theta,phi',[(0.,0.),(math.pi,1.),(.7,1.3),(math.pi/2,0.),(math.pi/2,math.pi/2)])
def test_basis_and_reconstruction_cartesian_reference(theta,phi):
    d=[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)]
    expected=basis_xyz(d)
    sh._basis_layout.cache_clear()
    first=sh._basis(3,theta,phi);second=sh._basis(3,theta,phi)
    assert first==second and first is not second
    assert first==pytest.approx(expected,abs=2e-15)
    coefficients=[[(-1)**i*(i+1)*.1,(i+1)*.2,0.] for i in range(9)]
    original=[row[:] for row in coefficients]
    reference=[math.fsum(row[c]*value for row,value in zip(coefficients,expected)) for c in range(3)]
    assert sh.reconstruct(coefficients,d)==pytest.approx(reference,abs=4e-15)
    assert coefficients==original


def grid(height=32,width=64):
    samples=[];weights=[];colors=[]
    for row in range(height):
        theta=math.pi*(row+.5)/height
        weight=(2*math.pi/width)*(math.cos(math.pi*row/height)-math.cos(math.pi*(row+1)/height))
        for col in range(width):
            phi=2*math.pi*(col+.5)/width
            d=[math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta)]
            sample=sh.SPHSample(theta,phi,Vector(3,d),9);sample.values[:]=basis_xyz(d)
            samples.append(sample);weights.append(weight);colors.append([2+d[0],3+d[1],4+d[2]])
    return samples,weights,colors


def test_constant_asymmetric_projection_and_diffuse_references():
    samples,weights,colors=grid()
    constant=sh.project_radiance(samples,[[2.,3.,4.]]*len(samples),weights)
    assert constant[0]==pytest.approx([math.sqrt(4*math.pi)*c for c in [2,3,4]],abs=2e-14)
    coefficients=sh.project_radiance(samples,colors,weights)
    a=math.sqrt(3/(4*math.pi))
    assert coefficients[1]==pytest.approx([0.,-1/a,0.],abs=.004)
    assert coefficients[2]==pytest.approx([0.,0.,1/a],abs=.004)
    assert coefficients[3]==pytest.approx([-1/a,0.,0.],abs=.004)
    original=[row[:] for row in coefficients]
    irradiance=sh.convolve_diffuse(coefficients)
    factors=[math.pi,2*math.pi/3,math.pi/4]
    assert irradiance==[[value*factors[math.isqrt(i)] for value in row] for i,row in enumerate(coefficients)]
    assert coefficients==original and all(a is not c for a,c in zip(irradiance,coefficients))
    for d in ([1.,0.,0.],[0.,1.,0.],[0.,0.,1.]):
        expected=[math.pi*ambient+2*math.pi/3*component for ambient,component in zip([2,3,4],d)]
        assert sh.reconstruct(irradiance,d)==pytest.approx(expected,abs=.013)


def test_compensated_projection_cancellation_and_channel_isolation():
    value=1/math.sqrt(4*math.pi)
    samples=[sh.SPHSample(0,0,Vector(3,[0,0,1]),1) for _ in range(3)]
    for sample in samples:sample.values[0]=value
    colors=[[float(2**54),0.,0.],[1.,0.,0.],[-float(2**54),0.,0.]]
    expected=float(sum(Fraction(row[0]*value) for row in colors))
    assert expected==value
    assert sh.project_radiance(samples,colors,[1.,1.,1.])==[[expected,0.,0.]]


def coefficients():
    # f(d)=[1+.2*x+.3*x*y, 2+.4*y*z, 3+.5*(x*x-y*y)]
    a=math.sqrt(3/(4*math.pi));c=math.sqrt(15/(4*math.pi));e=math.sqrt(15/(16*math.pi))
    rows=[[0.,0.,0.] for _ in range(9)]
    rows[0]=[math.sqrt(4*math.pi)*x for x in [1,2,3]]
    rows[3][0]=-.2/a;rows[4][0]=.3/c;rows[5][1]=-.4/c;rows[8][2]=.5/e
    return rows


def field(d):
    x,y,z=d
    return [1+.2*x+.3*x*y,2+.4*y*z,3+.5*(x*x-y*y)]


@pytest.mark.parametrize('direction',[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[.6,0.,.8]])
def test_active_rotation_and_noncommuting_composition(direction):
    rows=coefficients();snapshot=[r[:] for r in rows]
    qz=Quaternion([math.sqrt(.5),0.,0.,math.sqrt(.5)])
    qx=Quaternion([math.sqrt(.5),math.sqrt(.5),0.,0.])
    rotated=sh.rotate_coefficients(rows,qz)
    x,y,z=direction
    assert sh.reconstruct(rotated,direction)==pytest.approx(field([y,-x,z]),abs=6e-15)
    sequential=sh.rotate_coefficients(rotated,qx)
    combined=sh.rotate_coefficients(rows,qx*qz)
    assert sh.reconstruct(combined,direction)==pytest.approx(field([z,-x,-y]),abs=6e-15)
    for a,c in zip(sequential,combined):assert a==pytest.approx(c,abs=6e-15)
    for channel in range(3):
        scalar=sh.rotate_coefficients([r[channel] for r in rows],qz)
        assert scalar==[r[channel] for r in rotated]
        for start,end in ((0,1),(1,4),(4,9)):
            assert math.fsum(r[channel]**2 for r in rotated[start:end])==pytest.approx(math.fsum(r[channel]**2 for r in rows[start:end]),rel=3e-15,abs=3e-15)
    assert rows==snapshot and all(a is not c for a,c in zip(rows,rotated))


def test_layout_cache_is_bounded_and_higher_orders_unchanged():
    sh._basis_layout.cache_clear()
    for bands in range(1,21):
        assert sh._basis(bands,.73,1.27)==[sh.SPH(l,m,.73,1.27) for l in range(bands) for m in range(-l,l+1)]
    assert sh._basis_layout.cache_info().currsize==16
