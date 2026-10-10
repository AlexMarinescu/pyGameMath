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


def check_artifact(name, expected, actual):
    """Keep artifacts exact except eight audited unit-scale libm measurements.

    The allowance is one binary64 ULP at 1 between manifests, with a separate
    two-ULP-at-1 closed-form accuracy check. It is not a gem math tolerance.
    """
    import copy
    from decimal import Decimal, localcontext
    import json
    if name != 'measurements.json':
        if expected != actual:
            raise AssertionError('artifact mismatch: '+name)
        return
    def decode(payload):
        def pairs(items):
            result = {}
            for key, value in items:
                if key in result: raise AssertionError('duplicate measurement key')
                result[key] = value
            return result
        def invalid(value): raise AssertionError('nonfinite measurement: '+value)
        report = json.loads(payload, object_pairs_hook=pairs, parse_constant=invalid)
        canonical = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n').encode('utf-8')
        if payload != canonical: raise AssertionError('measurement serialization mismatch')
        check_numerics(report['scenes'])
        return report
    reference, observed = decode(expected), decode(actual)
    # These are exactly the fields affected by sin(pi/4)'s one-ULP platform
    # variation in the fixed 120-degree showcase. All other leaves stay exact.
    with localcontext() as context:
        context.prec = 90
        d = Decimal
        sqrt2, sqrt3, sqrt6 = d(2).sqrt(), d(3).sqrt(), d(6).sqrt()
        allowed = {
            (1, 'quaternion', 0): (sqrt6+sqrt2)/4,
            (1, 'axes', 0, 0): sqrt3/2,
            (1, 'axes', 0, 1): d('.5'),
            (1, 'axes', 1, 0): d('-.5'),
            (1, 'axes', 1, 1): sqrt3/2,
            (3, 'quaternion', 3): sqrt2/2,
            (3, 'axes', 0, 0): d(0),
            (3, 'axes', 1, 1): d(0),
        }
        difference_bound = math.ulp(1.0)
        accuracy_bound = d.from_float(2*math.ulp(1.0))
        normalized = copy.deepcopy(observed)
        for path, exact in allowed.items():
            a = reference['scenes']['quaternions']['frames']
            b = observed['scenes']['quaternions']['frames']
            target = normalized['scenes']['quaternions']['frames']
            for key in path[:-1]: a, b, target = a[key], b[key], target[key]
            av, bv = a[path[-1]], b[path[-1]]
            if type(av) is not float or type(bv) is not float or not all(map(math.isfinite, (av,bv))):
                raise AssertionError('invalid quaternion measurement')
            if abs(av-bv) > difference_bound:
                raise AssertionError('quaternion manifest difference exceeds audited bound')
            if any(abs(d.from_float(value)-exact) > accuracy_bound for value in (av,bv)):
                raise AssertionError('quaternion measurement fails closed-form reference')
            target[path[-1]] = av
        # Compare canonical bytes, keeping types, signed zeros, lengths, keys,
        # metadata, hashes and all unlisted numerical values exact.
        normalized_bytes = (json.dumps(normalized, indent=2, sort_keys=True, allow_nan=False)+'\n').encode('utf-8')
        if normalized_bytes != expected:
            raise AssertionError('unaudited measurement difference')
