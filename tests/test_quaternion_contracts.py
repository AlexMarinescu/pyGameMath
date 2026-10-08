"""Known answers distinguishing the legacy helper from rotation construction."""
import math
import pytest
from gem import quaternion as q, vector


@pytest.mark.parametrize('wrapped', [False, True])
@pytest.mark.parametrize('axis,angle,legacy,rotation', [
    ([0,0,2], 0, [0,0,0,1], [1,0,0,0]),
    ([0,0,2], 90, [0,0,0,1], [math.sqrt(.5),0,0,math.sqrt(.5)]),
    ([0,-3,0], -90, [0,0,-1,0], [math.sqrt(.5),0,math.sqrt(.5),0]),
    ([1,2,2], 120, [0,1/3,2/3,2/3], [.5,math.sqrt(3)/6,math.sqrt(3)/3,math.sqrt(3)/3]),
])
def test_independent_helper_results_and_ownership(wrapped,axis,angle,legacy,rotation):
    storage=axis[:]
    argument=vector.Vector(3,storage) if wrapped else storage
    old=q.quat_rotate_from_axis_angle(argument,angle)
    constructor=q.quat_from_axis_angle(argument,angle)
    assert old.data == pytest.approx(legacy,abs=1e-14)
    assert constructor.data == pytest.approx(rotation,abs=1e-14)
    assert isinstance(old,q.Quaternion) and isinstance(constructor,q.Quaternion)
    assert old is not constructor and old.data is not constructor.data
    assert old.data is not storage and constructor.data is not storage
    old.data[1]=12
    constructor.data[1]=13
    assert (argument.vector if wrapped else argument) is storage
    assert storage == axis


def test_documented_rotation_example():
    axis=[0,0,2]
    legacy=q.quat_rotate_from_axis_angle(axis,90)
    rotation=q.quat_from_axis_angle(axis,90)
    point=vector.Vector(3,[1,0,0]);storage=point.vector
    assert legacy.data == pytest.approx([0,0,0,1])
    assert q.quat_rotate_vector(rotation,point).vector == pytest.approx([0,1,0],abs=1e-14)
    assert point.vector is storage and storage == [1,0,0]
    assert axis == [0,0,2]


def test_constructor_retains_explicit_component_storage():
    data=[1,0,0,0]
    quat=q.Quaternion(data)
    assert quat.data is data
    copied=quat.conjugate()
    assert copied.data is not data
    data[1]=2
    assert quat.data == [1,2,0,0]
    assert copied.data == [1,0,0,0]
