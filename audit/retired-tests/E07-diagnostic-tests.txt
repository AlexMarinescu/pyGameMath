"""Diagnostic reproductions of unfinished E07 code, not a transport contract.

Test-only bypasses expose faults masked by the builtin-object typo. Production
code and its original strict expected failure remain unchanged.
"""
import math

import pytest

from gem.vector import Vector
from gem.experimental import sph_object, sph_sample


def fixture_geometry(count=1, allocate=False):
    vertices = [sph_object.SPHVertex(Vector(3, [i, 0, 0]),
                                   Vector(3, [0, 0, 1])) for i in range(count)]
    if allocate:
        for vertex in vertices:
            vertex.unshadowedCoeffs = [0.0]
            vertex.shadowedCoeffs = [0.0]
    obj = sph_object.SPHObject(list(range(count)), vertices)
    sample = sph_sample.SPHSample(0, 0, Vector(3, [0, 0, 1]), 1)
    sample.values[0] = 1 / math.sqrt(4 * math.pi)
    return obj, sample


def bypass_typo(monkeypatch, obj):
    monkeypatch.setattr(sph_object, 'object', [obj], raising=False)


def test_audit_original_builtin_subscript_failure():
    obj, sample = fixture_geometry()
    with pytest.raises(TypeError, match='not subscriptable'):
        sph_object.GenereateCoeffs(1, 1, [sample], [obj])


def test_audit_unallocated_arrays_after_typo_bypass(monkeypatch):
    obj, sample = fixture_geometry()
    bypass_typo(monkeypatch, obj)
    with pytest.raises(TypeError, match='item assignment'):
        sph_object.GenereateCoeffs(1, 1, [sample], [obj])


def test_audit_only_last_vertex_scaled_and_shadow_bias(monkeypatch):
    obj, sample = fixture_geometry(2, allocate=True)
    bypass_typo(monkeypatch, obj)
    assert sph_object.GenereateCoeffs(1, 1, [sample], [obj]) is None
    first, last = obj.vertices
    # Equal normals should have equal transfer; current output does not.
    assert first.unshadowedCoeffs == pytest.approx([1 / math.sqrt(4 * math.pi)])
    assert last.unshadowedCoeffs == pytest.approx([math.sqrt(4 * math.pi)])
    assert first.shadowedCoeffs == [0]
    assert last.shadowedCoeffs == pytest.approx([4 * math.pi])


@pytest.mark.parametrize('count, error', [(0, ZeroDivisionError), (2, IndexError)])
def test_audit_sample_count_errors(monkeypatch, count, error):
    obj, sample = fixture_geometry(allocate=True)
    bypass_typo(monkeypatch, obj)
    with pytest.raises(error):
        sph_object.GenereateCoeffs(count, 1, [sample], [obj])


def test_audit_empty_vertex_buffer(monkeypatch):
    obj, sample = fixture_geometry(0, allocate=True)
    bypass_typo(monkeypatch, obj)
    with pytest.raises(UnboundLocalError):
        sph_object.GenereateCoeffs(1, 1, [sample], [obj])


def test_audit_no_geometry_visibility_or_normal_validation(monkeypatch):
    obj, sample = fixture_geometry(allocate=True)
    bypass_typo(monkeypatch, obj)
    sph_object.GenereateCoeffs(1, 1, [sample], [obj])
    reference = obj.vertices[0].unshadowedCoeffs[:]
    obj.indices[:] = [999, -1]
    obj.vertices[0].position.vector[:] = [200, -300, 400]
    sph_object.GenereateCoeffs(1, 1, [sample], [obj])
    assert obj.vertices[0].unshadowedCoeffs == reference
    obj.vertices[0].normal.vector[2] = 2
    sph_object.GenereateCoeffs(1, 1, [sample], [obj])
    assert obj.vertices[0].unshadowedCoeffs == pytest.approx([2 * reference[0]])


def test_audit_constructor_aliases():
    position, normal = Vector(3, [1, 2, 3]), Vector(3, [0, 0, 1])
    vertex = sph_object.SPHVertex(position, normal)
    indices, vertices = [0], [vertex]
    obj = sph_object.SPHObject(indices, vertices)
    assert vertex.position is position and vertex.normal is normal
    assert obj.indices is indices and obj.vertices is vertices
    assert vertex.unshadowedCoeffs is None and vertex.shadowedCoeffs is None


def test_audit_empty_scene_noop():
    assert sph_object.GenereateCoeffs(0, 1, [], []) is None
