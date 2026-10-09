"""Independent coordinate references and structural SVG validation."""
import math
import xml.etree.ElementTree as ET


def close(actual, expected, tolerance=1e-12):
    if len(actual) != len(expected) or any(abs(a-b) > tolerance for a, b in zip(actual, expected)):
        raise AssertionError((actual, expected))


def segment_distance(p, a, b):
    dx, dy = b[0]-a[0], b[1]-a[1]
    denominator = dx*dx+dy*dy
    time = max(0., min(1., ((p[0]-a[0])*dx+(p[1]-a[1])*dy)/denominator)) if denominator else 0.
    return math.hypot(p[0]-a[0]-time*dx, p[1]-a[1]-time*dy)


def check_numerics(data):
    v = data['vectors']
    close(v['sum'], [2, 3, 0]); close(v['cross'], [0, 0, 7])
    close(v['normalized'], [3/math.sqrt(10), 1/math.sqrt(10), 0])
    assert v['dot'] == -1
    reference = [lambda x,y:(x,y), lambda x,y:(2*x,y), lambda x,y:(-y,x),
                 lambda x,y:(x+3,y-1), lambda x,y:(3-y,2*x-1), lambda x,y:(1-y,2*x+3)]
    for outputs, formula in zip(data['transforms']['results'], reference):
        for point, result in zip(data['transforms']['inputs'], outputs):
            close(result, [*formula(*point[:2]), 0, 1])
    for frame in data['quaternions']['frames']:
        angle = math.radians(120*frame['t'])
        close(frame['axes'][0], [math.cos(angle), math.sin(angle), 0])
        close(frame['axes'][1], [-math.sin(angle), math.cos(angle), 0])
        close(frame['quaternion'], [math.cos(angle/2), 0, 0, math.sin(angle/2)])
        assert abs(math.hypot(*frame['quaternion'])-1) < 1e-12
    b = data['bezier']
    close(b['quadratic_midpoint'], [2, 2]); close(b['cubic_midpoint'], [2, .875])
    errors = []
    for points, squared in zip(b['adaptive'], b['squared_tolerances']):
        close(points[0], [0, 0]); close(points[-1], [4, 1])
        assert all(a[0] < z[0] for a,z in zip(points, points[1:]))
        # Expanded cubic polynomials, independent of gem's Bernstein evaluator.
        maximum = 0.
        for i in range(1001):
            t = i/1000
            p = [3*t+3*t*t-2*t*t*t, 12*t-30*t*t+19*t*t*t]
            maximum = max(maximum, min(segment_distance(p,a,z) for a,z in zip(points,points[1:])))
        assert maximum <= math.sqrt(squared)+1e-12
        errors.append(maximum)
    assert len(b['adaptive'][1]) > len(b['adaptive'][0]) and errors[1] < errors[0]
    return dict(checks='vector, transform, SLERP, midpoint, endpoints, order and dense-reference distances',
                cubic_dense_reference_max_distances=errors)


def check_svg(text):
    root = ET.fromstring(text)
    assert root.tag == '{http://www.w3.org/2000/svg}svg'
    assert root.attrib['width'] == '960' and root.attrib['height'] == '560'
    assert root.find('{http://www.w3.org/2000/svg}title') is not None
    assert not any(node.tag.endswith(('script', 'image', 'foreignObject')) for node in root.iter())
    return [960, 560]
