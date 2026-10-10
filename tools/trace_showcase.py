"""Trace platform math used by the fixed quaternion showcase (diagnostics only)."""
import argparse
from decimal import Decimal, localcontext
import json
import math
from pathlib import Path
import platform
import struct
import sys

ROOT = Path(__file__).resolve().parents[1]
# Independent 100-digit mathematical constant, not binary64 math.pi.
PI = Decimal('3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679')


def sine(x):
    term = total = x
    for n in range(1, 200):
        term *= -x*x / Decimal((2*n)*(2*n+1))
        previous = total
        total += term
        if total == previous:
            return total
    raise AssertionError('reference sine did not converge')


def arctan(x):
    if x < 0:
        return -arctan(-x)
    if x > 1:
        return PI/2-arctan(1/x)
    if x > Decimal('.5'):
        return PI/4+arctan((x-1)/(x+1))
    term = total = x
    for n in range(1, 200):
        term *= -x*x
        previous = total
        total += term/Decimal(2*n+1)
        if previous == total:
            return total
    raise AssertionError('reference atan did not converge')


def reference(name, values):
    with localcontext() as context:
        context.prec = 90
        args = [Decimal.from_float(value) for value in values]
        if name == 'sin': return sine(args[0])
        if name == 'cos': return sine(PI/2-args[0])
        if name == 'hypot': return sum(x*x for x in args).sqrt()
        if name == 'radians': return args[0]*PI/180
        if name == 'atan2':
            y, x = args
            assert x > 0  # Every atan2 in this fixed showcase has positive X.
            return arctan(y/x)
        raise AssertionError(name)


def bits(value):
    return struct.unpack('>Q', struct.pack('>d', value))[0]


def describe(value, expected=None):
    result = {'value': value, 'hex': value.hex()}
    if expected is not None:
        rounded = float(expected)
        result.update(reference_decimal=str(expected), reference_hex=rounded.hex(),
                      ulps_from_rounded_reference=abs(bits(value)-bits(rounded)))
    return result


def trace():
    sys.path.insert(0, str(ROOT))
    from gem.quaternion import Quaternion, quat_from_axis_angle, quat_slerp, quat_rotate_vector
    from gem.vector import Vector
    calls = []
    names = ('radians', 'sin', 'cos', 'hypot', 'atan2')
    originals = {name: getattr(math, name) for name in names}
    def wrapper(name):
        def call(*args):
            value = originals[name](*args)
            calls.append({'operation': name, 'arguments': [describe(float(a)) for a in args],
                          'result': describe(value, reference(name, list(map(float, args))))})
            return value
        return call
    frames = []
    try:
        for name in names: setattr(math, name, wrapper(name))
        q0, q1 = Quaternion(), quat_from_axis_angle([0., 0., 1.], 120)
        endpoint = [describe(x) for x in q1.data]
        for t in (0., .25, .5, .75, 1.):
            start = len(calls)
            q = quat_slerp(q0, q1, t)
            axes = [quat_rotate_vector(q, Vector(3, list(axis))).vector
                    for axis in ((1., 0., 0.), (0., 1., 0.))]
            with localcontext() as context:
                context.prec = 90
                angle = PI*Decimal.from_float(t)*Decimal(2)/3
                s, c = sine(angle), sine(PI/2-angle)
                expected = ((c, s, Decimal(0)), (-s, c, Decimal(0)))
                frames.append({'t': t, 'math_call_range': [start, len(calls)],
                               'quaternion': [describe(x) for x in q.data],
                               'axes': [[describe(x, e) for x, e in zip(axis, ref)]
                                        for axis, ref in zip(axes, expected)]})
    finally:
        for name, function in originals.items(): setattr(math, name, function)
    return {'schema_version': 1, 'python': sys.version, 'implementation': platform.python_implementation(),
            'compiler': platform.python_compiler(), 'platform': platform.platform(),
            'machine': platform.machine(), 'libc': platform.libc_ver(),
            'float_info': str(sys.float_info), 'reference_precision': 90,
            'endpoint': endpoint, 'math_calls': calls, 'frames': frames}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(trace(), indent=2, allow_nan=False)+'\n', encoding='utf-8')


if __name__ == '__main__': main()
