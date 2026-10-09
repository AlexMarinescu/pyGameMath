"""Independent numerical and dispatch checks for focused allocation reductions."""
import math
import struct
from decimal import Decimal, localcontext
from fractions import Fraction
import pytest
from gem import vector as v, quaternion as q


def axis(degrees):
    angle=math.radians(degrees)/2
    return q.Quaternion([math.cos(angle),0.,0.,math.sin(angle)])


def hamilton(a,b):
    # Scalar/vector definition, independent of the component kernel.
    aw,*av=map(Fraction,a);bw,*bv=map(Fraction,b)
    cross=[av[1]*bv[2]-av[2]*bv[1],av[2]*bv[0]-av[0]*bv[2],av[0]*bv[1]-av[1]*bv[0]]
    return [aw*bw-sum(x*y for x,y in zip(av,bv))]+[aw*y+bw*x+c for x,y,c in zip(av,bv,cross)]


@pytest.mark.parametrize('size',[0,1,2,3,4,8])
@pytest.mark.parametrize('scale',[0.,1.,1e-300,1e300,math.ldexp(1.,-1074)])
def test_normalization_decimal_reference_and_ownership(size,scale):
    values=[(-1.)**i*(i+1)*scale for i in range(size)]
    obj=v.Vector(size,values)
    with localcontext() as context:
        context.prec=800
        data=[Decimal.from_float(x) for x in values]
        length=sum((x*x for x in data),Decimal(0)).sqrt()
        expected=[float(x/length) for x in data] if length else [0.]*size
    result=obj.normalize()
    assert result.vector==pytest.approx(expected,rel=2e-15,abs=math.ldexp(1.,-1074))
    assert result.vector is not values and obj.vector is values
    snapshot=values[:]
    assert obj.i_normalize() is obj
    assert obj.vector==result.vector and obj.vector is not values and values==snapshot
    if size==4:
        quat=q.Quaternion(values)
        result=quat.normalize()
        assert result.data==pytest.approx(expected if length else [1,0,0,0],rel=2e-15,abs=math.ldexp(1.,-1074))
        assert result.data is not values and quat.data is values
        assert quat.i_normalize() is quat and quat.data==result.data and values==snapshot


@pytest.mark.parametrize('data',[[1e308,-1e308,1e308,1e308],[1e300,1e-300,-1.,0.],
                                  [1.,1e-300,-1e-300,0.],[-0.,0.,-0.,0.]])
def test_overrange_and_mixed_normalization(data):
    with localcontext() as context:
        context.prec=800
        exact=[Decimal.from_float(x) for x in data]
        norm=sum(x*x for x in exact).sqrt()
        expected=[float(x/norm) for x in exact] if norm else [1,0,0,0]
    assert q.Quaternion(data).normalize().data==pytest.approx(expected,rel=2e-15,abs=math.ldexp(1.,-1074))


@pytest.mark.parametrize('values',[[2.,-1.,3.,.5],[.5,.5,.5,.5],[-1.,0.,0.,0.]])
@pytest.mark.parametrize('point',[[1.,2.,3.],[-2.,.5,4.],[0.,0.,0.]])
def test_rotation_fraction_sandwich_and_nonunit_inverse(values,point):
    quat=q.Quaternion(values);vec=v.Vector(3,point)
    conjugate=[values[0]]+[-x for x in values[1:]]
    reference=hamilton(hamilton(values,[0]+point),conjugate)
    result=q.quat_rotate_vector(quat,vec)
    assert result.vector==pytest.approx(list(map(float,reference[1:])),rel=2e-15,abs=1e-15)
    assert result.vector is not point and quat.data is values and vec.vector is point
    norm2=sum(Fraction(x)**2 for x in values)
    expected=[float(Fraction(x)/norm2) for x in conjugate]
    assert quat.inverse().data==pytest.approx(expected,rel=2e-15,abs=1e-15)
    assert (quat*quat.inverse()).data==pytest.approx([1,0,0,0],abs=1e-15)
    assert (quat.inverse()*quat).data==pytest.approx([1,0,0,0],abs=1e-15)


@pytest.mark.parametrize('t',[0.,.13,.5,.87,1.])
@pytest.mark.parametrize('angle',[1e-12,.001,90.,179.])
@pytest.mark.parametrize('sign',[1.,-1.])
def test_shortest_slerp_independent_angles(angle,t,sign):
    first=axis(0);last=axis(angle);last.data=[sign*x for x in last.data]
    original=(first.data[:],last.data[:])
    result=first.slerp(last,t)
    assert result.data==pytest.approx(axis(angle*t).data,abs=1e-12,rel=0)
    assert math.hypot(*result.data)==pytest.approx(1.,abs=1e-12,rel=0)
    assert (first.data,last.data)==original
    assert result.data is not first.data and result.data is not last.data


@pytest.mark.parametrize('t',[0.,.2,.5,.8,1.])
def test_squad_references_and_ownership(t):
    controls=[axis(x) for x in (0,180,90)]
    storage=[obj.data for obj in controls];snapshots=[data[:] for data in storage]
    legacy=q.quat_squad(*controls,t)
    blend=2*t*(1-t)
    # The historical final blend is linear when cos(45*t degrees) >= .95.
    if math.cos(math.radians(45*t)) >= .95:
        expected=[(1-blend)*x+blend*y for x,y in zip(axis(90*t).data,axis(180*t).data)]
    else:
        expected=axis(90*t+180*t*t*(1-t)).data
    assert legacy.data==pytest.approx(expected,abs=1e-12)
    standard=[axis(x) for x in (0,90,30,120)]
    result=q.squad4(*standard,t)
    assert result.data==pytest.approx(axis(90*t+60*t*(1-t)).data,abs=1e-12)
    assert math.hypot(*result.data)==pytest.approx(1,abs=1e-12)
    assert [obj.data for obj in controls]==snapshots
    assert all(obj.data is data for obj,data in zip(controls,storage))
    assert all(result.data is not obj.data for obj in standard)


def test_cross_known_answer_and_no_alias():
    a=[2.,-3.,4.];b=[-5.,6.,7.]
    result=v.cross(v.Vector(3,a),v.Vector(3,b))
    assert result.vector==[-45.,-34.,-3.]
    assert a==[2,-3,4] and b==[-5,6,7] and result.vector is not a and result.vector is not b


def test_conjugation_signs_and_in_place_storage():
    data=[-0.,0.,-0.,2.]
    quat=q.Quaternion(data);result=quat.conjugate()
    expected=[-0.,-0.,0.,-2.]
    assert [struct.pack('d',x) for x in result.data]==[struct.pack('d',x) for x in expected]
    assert result.data is not data
    assert quat.i_conjugate() is quat and quat.data==result.data and quat.data is not data
    assert data==[-0.,0.,-0.,2.]


def test_subclass_operator_dispatch_is_retained():
    events=[]
    class Custom(q.Quaternion):
        def __mul__(self,other):
            events.append('multiply')
            return super().__mul__(other)
        def conjugate(self):
            events.append('conjugate')
            return super().conjugate()
    obj=Custom([1.,0.,0.,0.]);other=Custom([0.,0.,0.,1.])
    assert q.quat_rotate_vector(obj,v.Vector(3,[1,2,3])).vector==[1,2,3]
    assert events==['multiply','conjugate']
    events.clear()
    assert q.quat_slerp(obj,other,.5).data==pytest.approx([math.sqrt(.5),0,0,math.sqrt(.5)])
    assert events==['multiply','multiply']


@pytest.mark.parametrize('values',[[float('nan'),0.,1.],[float('inf'),0.,1.],[-float('inf'),1.,2.]])
def test_nonfinite_normalization_preserves_arithmetic(values):
    out=v.normalize(3,values)
    if math.isnan(values[0]):assert all(math.isnan(x) for x in out)
    else:assert math.isnan(out[0]) and out[1:]==[0.,0.]


@pytest.mark.parametrize('operation',[
    lambda:q.quat_lerp(axis(0),axis(90),0),
    lambda:q.quat_rotate_vector(axis(0),[1,2,3]),
    lambda:q.quat_rotate_vector(axis(0),v.Vector(2,[1,2])),
    lambda:q.quat_conjugate([1.,2.,3.]),
    lambda:v.normalize(3,[1.,2.]),
    lambda:q.Quaternion([0.,0.,0.,0.]).inverse(),
])
def test_existing_invalid_domains(operation):
    with pytest.raises((TypeError,IndexError,ZeroDivisionError)):
        operation()


def test_transform_and_clamp_contracts_remain_intact():
    values=[1.,2.,3.];rows=[[0.,1.,0.,0.],[-1.,0.,0.,0.],[0.,0.,1.,0.],[4.,-5.,6.,1.]]
    obj=v.Vector(3,values)
    result=obj.transform(values,rows)
    assert result.vector==[2.,-4.,9.] and values==[1.,2.,3.]
    assert obj.i_transform(values,rows) is obj and obj.vector==result.vector and values==[1,2,3]
    lo=[0.,0.,0.];hi=[2.,2.,2.]
    result=v.clamp(3,values,lo,hi)
    assert result.vector==[1,2,2] and result.vector is not values
    assert obj.i_clamp(3,values,lo,hi) is obj and obj.vector==[1,2,2] and values==[1,2,3]
