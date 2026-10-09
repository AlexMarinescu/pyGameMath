"""Graphics geometry comes from supported gem APIs; SVG is presentation only."""
from gem.vector import Vector, cross
from gem.matrix import Matrix
from gem.quaternion import Quaternion, quat_from_axis_angle, quat_rotate_vector, quat_slerp
from gem.bezier import quadraticBezierPoint, cubicBezierPoint, BezierPath
from .svg import Canvas, chart, BLUE, ORANGE, GREEN, MUTED


def vectors():
    a, b = Vector(3, [3., 1., 0.]), Vector(3, [-1., 2., 0.])
    total, unit = a+b, a.normalize()
    data = dict(a=a.vector, b=b.vector, sum=total.vector, normalized=unit.vector,
                dot=a.dot(b), cross=cross(a, b).vector)
    c = Canvas('Vector geometry', 'Displacement, length and orientation from the same pair of vectors')
    screen = chart(c, (185, 405), 78, (-1.5, 4, -1, 3.7))
    origin = screen((0, 0))
    c.arrow(origin, screen(a.vector), BLUE)
    c.arrow(origin, screen(b.vector), ORANGE)
    c.arrow(screen(a.vector), screen(total.vector), ORANGE)
    c.arrow(origin, screen(total.vector), GREEN)
    c.text(425, 336, 'a = (3, 1)', color=BLUE)
    c.text(42, 218, 'b = (-1, 2)', color=ORANGE)
    c.text(365, 150, 'a + b = (2, 3)', color=GREEN)
    c.text(423, 210, 'displaced b', 14, ORANGE)
    c.text(600, 145, 'Normalize without changing a', 19)
    c.arrow((620, 265), (620+150*unit.vector[0], 265-150*unit.vector[1]), BLUE)
    c.text(600, 305, 'unit a = (0.948683, 0.316228, 0)', 15)
    c.text(600, 365, 'a dot b = -1  (obtuse angle)', 18)
    c.text(600, 405, 'a cross b = (0, 0, 7)', 18, GREEN)
    c.text(600, 433, '+Z points out of the XY diagram', 15)
    c.text(32, 530, 'Arrow positions are coordinates; normalized direction has unit length.', 15, MUTED)
    return c.finish(), data


def transforms():
    points = [Vector(4, [x, y, 0., 1.]) for x, y in [(0, 0), (1, 0), (1, 2), (0, 2)]]
    s = Matrix(4).scale(Vector(3, [2., 1., 1.]))
    r = Matrix(4).rotate(Vector(3, [0., 0., 1.]), 90)
    t = Matrix(4).translate(Vector(3, [3., -1., 0.]))
    models = [Matrix(4), s, r, t, s*r*t, s*t*r]
    names = ['Original', 'Scale X by 2', 'Rotate +90 degrees', 'Translate (3, -1)', 'S then R then T', 'S then T then R']
    results = [[(m*p).vector for p in points] for m in models]
    c = Canvas('Transformation order', 'Matrix * Vector evaluates the row product vM; A * B applies A, then B')
    for i, (name, result) in enumerate(zip(names, results)):
        column, row = i % 3, i // 3
        left, top = 32+column*308, 112+row*205
        c.text(left, top, name, 17)
        screen = chart(c, (left+88, top+145), 23, (-2, 5, -2, 4))
        c.path([screen(p.vector) for p in points], '#a4b3bf', 2, True, True)
        c.path([screen(p) for p in result], BLUE if i < 4 else ORANGE, 3, True)
        for p in result:
            c.circle(screen(p), 3, BLUE if i < 4 else ORANGE)
    c.text(32, 535, 'Combined: (x, y) -> (3 - y, 2x - 1). Reversing T and R rotates the translation.', 15, MUTED)
    return c.finish(), dict(inputs=[p.vector for p in points], labels=names, results=results)


def quaternions():
    q0, q1 = Quaternion(), quat_from_axis_angle([0., 0., 1.], 120)
    c = Canvas('Quaternion orientation with SLERP', 'Unit orientations; active +Z rotation moves +X toward +Y')
    data = []
    for i, time in enumerate([0., .25, .5, .75, 1.]):
        q = quat_slerp(q0, q1, time)
        axes = [quat_rotate_vector(q, Vector(3, list(axis))).vector for axis in [(1., 0., 0.), (0., 1., 0.)]]
        origin = (108+i*184, 300)
        c.circle(origin, 4, MUTED)
        for axis, color in zip(axes, [BLUE, ORANGE]):
            c.arrow(origin, (origin[0]+64*axis[0], origin[1]-64*axis[1]), color)
        c.text(origin[0]-45, 410, 't = '+str(time), 17)
        c.text(origin[0]-45, 438, str(round(time*120))+' degrees', 16)
        data.append(dict(t=time, quaternion=q.data, axes=axes))
    c.text(32, 155, 'Blue: local X     Orange: local Y     All frames share world +Z', 18)
    c.text(32, 515, 'SLERP follows the shortest orientation path; samples are uniform in t, not animation keyframes.', 15, MUTED)
    return c.finish(), dict(endpoint=q1.data, frames=data)


def bezier():
    quadratic = [Vector(2, list(p)) for p in [(0., 0.), (2., 4.), (4., 0.)]]
    cubic = [Vector(2, list(p)) for p in [(0., 0.), (1., 4.), (3., -2.), (4., 1.)]]
    qpoints = [quadraticBezierPoint(i/64, *quadratic).vector for i in range(65)]
    cpoints = [cubicBezierPoint(i/128, *cubic).vector for i in range(129)]
    path = BezierPath()
    path.setControlPoints(cubic)
    samples = []
    for tolerance in [.16, .0025]:
        path.minimum_sqr_distance = tolerance
        samples.append([p.vector for p in path.findDrawingPoints(0)])
    c = Canvas('Bezier evaluation and adaptive sampling', 'Control polygons are dashed; uniform parameter steps are not distance-uniform traversal')
    for i, (title, controls, curve) in enumerate([('Quadratic', quadratic, qpoints), ('Cubic', cubic, cpoints)]):
        left = 42+470*i
        screen = chart(c, (left+65, 345), 55, (-.2, 4.5, -1.2, 3.4))
        c.text(left, 112, title, 20)
        c.path([screen(p.vector) for p in controls], MUTED, 2, dashed=True)
        c.path([screen(p) for p in curve], BLUE, 3)
        for index, p in enumerate(controls):
            c.circle(screen(p.vector), 5, MUTED)
            x, y = screen(p.vector)
            # Place labels below the upper control, avoiding the panel heading.
            c.text(x+8, max(142, y+(42 if index == len(controls)-1 else 20)), 'P'+str(index), 14)
        if i == 0:
            for j in range(9):
                c.circle(screen(quadraticBezierPoint(j/8, *quadratic).vector), 4, ORANGE)
            c.text(left, 486, 'Orange: uniform t = 0, 1/8, ..., 1', 15, ORANGE)
        else:
            for p in samples[0]:
                c.circle(screen(p), 6, ORANGE)
            for p in samples[1]:
                c.circle(screen(p), 2.5, GREEN)
            c.text(left, 486, 'distance 0.4: '+str(len(samples[0]))+' points (orange)', 15, ORANGE)
            c.text(left, 511, 'distance 0.05: '+str(len(samples[1]))+' points (green)', 15, GREEN)
    c.text(42, 540, 'BezierPath uses squared tolerances 0.16 and 0.0025; subdivision stops at depth 16.', 14, MUTED)
    return c.finish(), dict(quadratic_controls=[p.vector for p in quadratic], cubic_controls=[p.vector for p in cubic],
                           quadratic_midpoint=quadraticBezierPoint(.5, *quadratic).vector,
                           cubic_midpoint=cubicBezierPoint(.5, *cubic).vector,
                           squared_tolerances=[.16, .0025], adaptive=samples)


SCENES = dict(vectors=vectors, transforms=transforms, quaternions=quaternions, bezier=bezier)
