"""Phase 2C regressions; known angles and independent unit checks."""
import math
import random
import pytest
from gem import common, vector


@pytest.mark.parametrize('degrees,radians', [
    (0, 0), (30, math.pi/6), (45, math.pi/4), (60, math.pi/3),
    (90, math.pi/2), (180, math.pi), (270, 3*math.pi/2), (360, 2*math.pi),
    (-30, -math.pi/6), (-90, -math.pi/2), (-180, -math.pi),
    (-360, -2*math.pi), (720, 4*math.pi),
])
def test_known_angles(degrees, radians):
    assert common.radiansToDegrees(radians) == pytest.approx(degrees, rel=1e-14, abs=1e-14)
    assert common.degreesToRadians(degrees) == pytest.approx(radians, rel=1e-14, abs=1e-14)


@pytest.mark.parametrize('value', [-1e6, -720, -1, -1e-8, 0, 1e-8, 1, 720, 1e6])
def test_conversion_roundtrips(value):
    assert common.degreesToRadians(common.radiansToDegrees(value)) == pytest.approx(value, rel=1e-14, abs=0)
    assert common.radiansToDegrees(common.degreesToRadians(value)) == pytest.approx(value, rel=1e-14, abs=0)


@pytest.mark.parametrize('seed', range(20))
def test_conversion_against_standard_library(seed):
    rng = random.Random(seed)
    degrees, radians = rng.uniform(-10000, 10000), rng.uniform(-100, 100)
    assert common.degreesToRadians(degrees) == pytest.approx(math.radians(degrees), rel=1e-14)
    assert common.radiansToDegrees(radians) == pytest.approx(math.degrees(radians), rel=1e-14)
    assert common.degreesToRadians(-degrees) == -common.degreesToRadians(degrees)
    assert common.radiansToDegrees(-radians) == -common.radiansToDegrees(radians)


def test_legacy_keyword_parameter_names():
    assert common.radiansToDegrees(degrees=math.pi) == pytest.approx(180)
    assert common.degreesToRadians(radians=180) == pytest.approx(math.pi)


def V(values):
    return vector.Vector(len(values), list(values))


@pytest.mark.parametrize('ratio', [0.5, 2/3, 1.0, 1.5, 2.0])
def test_refraction_normal_incidence(ratio):
    incident, normal = V([0,-1,0]), V([0,1,0])
    result = vector.refract(ratio, incident, normal)
    assert result.vector == pytest.approx([0,-1,0], abs=1e-14)
    assert result is not incident and result is not normal
    assert incident.vector == [0,-1,0] and normal.vector == [0,1,0]


@pytest.mark.parametrize('ratio,incident,normal,expected', [
    (0.5, [0.6,-0.8,0], [0,1,0], [0.3,-math.sqrt(0.91),0]),
    (1.5, [0.6,-0.8,0], [0,1,0], [0.9,-math.sqrt(0.19),0]),
    (1.0, [0.6,-0.8,0], [0,1,0], [0.6,-0.8,0]),
    (1.5, [1,0,0], [-0.8,0.6,0],
     [0.54+0.8*math.sqrt(0.19), 0.72-0.6*math.sqrt(0.19), 0]),
])
def test_refraction_oblique_known_answers(ratio, incident, normal, expected):
    incoming, surface = V(incident), V(normal)
    result = vector.refract(IOR=ratio, incidentVec=incoming, normal=surface)
    assert isinstance(result, vector.Vector)
    assert result.size == 3
    assert result.vector == pytest.approx(expected, abs=1e-14)
    assert incoming.vector == incident and surface.vector == normal


@pytest.mark.parametrize('ratio', [1.25, 1.5, 2.0])
def test_refraction_total_internal_reflection(ratio):
    # Pick an angle beyond asin(n2/n1), without using the kernel's k formula.
    theta = (math.asin(1/ratio) + math.pi/2) / 2
    incident, normal = V([math.sin(theta),-math.cos(theta),0]), V([0,1,0])
    original = list(incident.vector)
    result = vector.refract(ratio, incident, normal)
    assert result.size == 3 and result.vector == [0,0,0]
    result.vector[0] = 99
    assert incident.vector == original and normal.vector == [0,1,0]


def test_refraction_exact_critical_boundary():
    # eta=1.25 and dot=-0.6 evaluate k to exactly zero in float arithmetic.
    incident, normal = V([0.8,-0.6,0]), V([0,1,0])
    result = vector.refract(1.25, incident, normal)
    assert result.vector == pytest.approx([1,0,0], abs=1e-14)


@pytest.mark.parametrize('n1,n2', [
    pytest.param(1.0,1.5,id='air-to-glass'),
    pytest.param(1.5,1.0,id='glass-to-air'),
])
@pytest.mark.parametrize('theta', [0, math.pi/6], ids=['normal','30-degrees'])
def test_refraction_between_air_and_glass(n1, n2, theta):
    incident = V([math.sin(theta),-math.cos(theta),0])
    transmitted_angle = math.asin((n1/n2)*math.sin(theta))
    expected = [math.sin(transmitted_angle),-math.cos(transmitted_angle),0]
    result = vector.refract(n1/n2, incident, V([0,1,0]))
    assert result.vector == pytest.approx(expected, abs=1e-14)
    assert n1*math.sin(theta) == pytest.approx(n2*result.vector[0], abs=1e-14)


@pytest.mark.parametrize('seed', range(20))
def test_refraction_snell_law_and_reversibility(seed):
    rng = random.Random(seed)
    n1, n2 = rng.uniform(1,2.5), rng.uniform(1,2.5)
    ratio = n1/n2
    # Independent geometry: a unit normal and perpendicular unit tangent.
    azimuth, polar = rng.uniform(0,2*math.pi), rng.uniform(0.2,math.pi-0.2)
    normal = [math.sin(polar)*math.cos(azimuth),
              math.sin(polar)*math.sin(azimuth), math.cos(polar)]
    tangent = [-math.sin(azimuth), math.cos(azimuth), 0]
    limit = math.asin(min(1,1/ratio))
    theta_i = rng.uniform(0.1,0.8) * limit
    incident = [math.sin(theta_i)*t-math.cos(theta_i)*n for t,n in zip(tangent,normal)]
    theta_t = math.asin(ratio*math.sin(theta_i))
    expected = [math.sin(theta_t)*t-math.cos(theta_t)*n for t,n in zip(tangent,normal)]
    incoming, surface = V(incident), V(normal)
    result = vector.refract(ratio, incoming, surface)
    assert result.vector == pytest.approx(expected, abs=1e-13)
    assert sum(x*x for x in result.vector) == pytest.approx(1, abs=1e-13)
    transmitted_sine = sum(x*t for x,t in zip(result.vector,tangent))
    assert n1*math.sin(theta_i) == pytest.approx(n2*transmitted_sine, abs=1e-13)
    assert sum(x*n for x,n in zip(result.vector,normal)) < 0
    reversed_ray = vector.refract(1/ratio, V([-x for x in result.vector]), V([-x for x in normal]))
    assert reversed_ray.vector == pytest.approx([-x for x in incident], abs=1e-13)
    assert incoming.vector == incident and surface.vector == normal
