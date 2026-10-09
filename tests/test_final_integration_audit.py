"""Independent cross-module checks: rational affine maps and analytical rotations.

No gem conversion, inverse, SH basis or interpolation supplies the sole oracle.
"""
from fractions import Fraction as F
import ctypes
import math
import random

import pytest

from gem import bezier, common, legendre, matrix, plane, quaternion, ray
from gem import spherical_harmonics as sh
from gem.vector import Vector, reflect, refract


def V(values):
    return Vector(len(values), list(values))


def unit(values):
    length = math.hypot(*values)
    return [x/length for x in values]


def rotate(values, axis, degrees):
    axis = unit(axis)
    c, s = math.cos(math.radians(degrees)), math.sin(math.radians(degrees))
    dot = math.fsum(x*y for x,y in zip(axis, values))
    cross = [axis[1]*values[2]-axis[2]*values[1],
             axis[2]*values[0]-axis[0]*values[2], axis[0]*values[1]-axis[1]*values[0]]
    return [c*v+s*w+(1-c)*dot*n for v,w,n in zip(values, cross, axis)]


def row(values, rows):
    return [sum((F(values[i])*F(rows[i][j]) for i in range(len(values))), F())
            for j in range(len(rows))]


def inverse_reference(rows):
    n = len(rows)
    augmented = [[F(x) for x in r]+[F(i == j) for j in range(n)] for i,r in enumerate(rows)]
    for column in range(n):
        pivot = next(i for i in range(column,n) if augmented[i][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [x/divisor for x in augmented[column]]
        for i in range(n):
            if i != column:
                factor = augmented[i][column]
                augmented[i] = [x-factor*y for x,y in zip(augmented[i],augmented[column])]
    return [r[n:] for r in augmented]


def bernstein(points, t):
    degree = len(points)-1
    weights = [math.comb(degree,i)*F(t)**i*(1-F(t))**(degree-i) for i in range(degree+1)]
    return [sum((F(p[j])*w for p,w in zip(points, weights)), F()) for j in range(len(points[0]))]


def basis(direction):
    x,y,z = direction
    a,b,c,d = math.sqrt(3/(4*math.pi)), math.sqrt(15/(4*math.pi)), math.sqrt(5/(16*math.pi)), math.sqrt(15/(16*math.pi))
    return [1/math.sqrt(4*math.pi),-a*y,a*z,-a*x,b*x*y,-b*y*z,c*(3*z*z-1),-b*x*z,d*(x*x-y*y)]


def evaluate(coefficients, direction):
    values = basis(direction)
    return [math.fsum(c[channel]*value for c,value in zip(coefficients,values)) for channel in range(3)]


@pytest.mark.parametrize('seed', range(12))
@pytest.mark.parametrize('scale', [1e-300, 1., 1e300])
def test_rotation_construction_conversion_and_composition(seed, scale):
    rng = random.Random(475000+seed)
    axis = [rng.uniform(-2,2) for _ in range(3)]
    values = [rng.uniform(-3,3) for _ in range(3)]
    degrees = [-180., -90., -.01, 0., .01, 90., 180., 73.][seed % 8]
    storage = [x*scale for x in axis]
    axis_vector = Vector(3, storage)
    q = quaternion.quat_from_axis_angle(axis_vector, degrees)
    rotation = matrix.Matrix(4).rotate(axis_vector, degrees)
    expected = rotate(values, axis, degrees)
    for actual in (quaternion.quat_rotate_vector(q,V(values)),
                   (rotation*V(values+[0.])).xyz(), (q.toMatrix()*V(values+[0.])).xyz(),
                   quaternion.quat_rotate_vector(quaternion.quat_from_matrix(rotation), V(values))):
        assert actual.vector == pytest.approx(expected, abs=2e-14, rel=2e-14)
    other_axis, other_angle = [2.,-1.,3.], -47.
    q2 = quaternion.quat_from_axis_angle(other_axis, other_angle)
    expected2 = rotate(expected, other_axis, other_angle)
    assert quaternion.quat_rotate_vector(q2*q,V(values)).vector == pytest.approx(expected2, abs=3e-14)
    composed = rotation*q2.toMatrix()
    assert (composed*V(values+[1.])).vector == pytest.approx(expected2+[1.], abs=3e-14)
    inv = composed.inverse()
    for i,r in enumerate(inv.matrix):
        assert r == pytest.approx([float(x) for x in inverse_reference(composed.matrix)[i]], abs=4e-14)
    exported = q.toMatrix()
    assert common.list_2d_to_1d([list(r) for r in exported.c_matrix]) == [
        ctypes.c_float(x).value for x in common.list_2d_to_1d(exported.matrix)]
    assert axis_vector.vector is storage and storage == [x*scale for x in axis]


@pytest.mark.parametrize('seed', range(12))
def test_rotated_reflection_refraction_and_ray_geometry(seed):
    axis, degrees = [1.,-2.,2.], 13.*seed-71.
    incident, normal = [.6,-.8,0.], [0.,1.,0.]
    rotated_i, rotated_n = rotate(incident,axis,degrees), rotate(normal,axis,degrees)
    q = quaternion.quat_from_axis_angle(axis,degrees)
    saved = q.data[:]
    for eta in (2./3.,1.,1.5):
        transmitted = [eta*.6, -math.sqrt(1.-(eta*.6)**2), 0.]
        assert refract(eta,V(rotated_i),V(rotated_n)).vector == pytest.approx(
            rotate(transmitted,axis,degrees), abs=3e-14)
    assert reflect(V(rotated_i),V(rotated_n)).vector == pytest.approx(rotate([.6,.8,0.],axis,degrees),abs=3e-14)
    start, direction = V([1.,2.,-3.]), V([3.,-4.,0.])
    original_direction_storage = direction.vector
    r = ray.Ray(start,direction)
    assert r.start is start and r.dir is direction and r.distance == 5.
    assert original_direction_storage == [3.,-4.,0.] and direction.vector is not original_direction_storage
    end = r.end
    copy = r.duplicate()
    assert r.rotateUsingQuaternion(q) is None
    assert r.translate(matrix.Matrix(4).translate(V([2.,-1.,4.]))) is None
    expected_start = [x+y for x,y in zip(rotate([1.,2.,-3.],axis,degrees),[2.,-1.,4.])]
    expected_dir = rotate(incident,axis,degrees)
    assert r.start.vector == pytest.approx(expected_start,abs=2e-14)
    assert r.dir.vector == pytest.approx(expected_dir,abs=2e-14)
    for t in (-2.,0.,5.):
        assert (r.start+r.dir*t).vector == pytest.approx([s+t*d for s,d in zip(expected_start,expected_dir)],abs=5e-14)
    assert r.end is end and r.distance == 5. and q.data == saved
    assert start.vector == [1.,2.,-3.] and copy.start.vector == start.vector
    assert copy.dir.vector == incident and copy.end is not end


@pytest.mark.parametrize('sx', [-2., .5, 2.])
@pytest.mark.parametrize('sz', [-.5, 1., 4.])
def test_affine_plane_normal_uses_inverse_transpose_and_winding(sx,sz):
    rows = [[sx,.5,0.,0.],[0.,2.,.25,0.],[0.,0.,sz,0.],[3.,-2.,1.,1.]]
    transform = matrix.Matrix(4, rows)
    inverse = inverse_reference(rows)
    original_normal = [-1.,1.,1.]
    # Covectors use inverse applied as a column, not the position transform.
    dual = [sum((inverse[i][j]*F(original_normal[j]) for j in range(3)),F()) for i in range(3)]
    expected_d = F(-2)-sum(F(rows[3][i])*dual[i] for i in range(3))
    normal_matrix = matrix.Matrix(3, common.convertM4to3(transform.inverse().matrix)).transpose()
    assert (normal_matrix*V(original_normal)).vector == pytest.approx([float(x) for x in dual],abs=2e-14)
    p = plane.Plane(); p.fromCoeffs(*map(float,dual),float(expected_d))
    controls = [[0.,0.,2.],[2.,0.,4.],[0.,2.,0.]]
    transformed = []
    for point in controls:
        expected = [float(x) for x in row(point+[1.],rows)]
        actual = Vector(3).transform(point,rows)
        assert actual.vector == expected[:3]
        assert p.dot(V(expected)) == pytest.approx(0.,abs=3e-14)
        transformed.append(actual)
    from_points = plane.Plane(); from_points.fromPoints(*transformed)
    winding = -1 if sx*sz < 0 else 1
    assert from_points.normal.vector == pytest.approx(unit([winding*float(x) for x in dual]),abs=2e-14)
    assert from_points.bestFitNormal(transformed).vector == pytest.approx(from_points.normal.vector,abs=2e-14)


@pytest.mark.parametrize('seed', range(8))
@pytest.mark.parametrize('perspective', [False,True])
def test_camera_pose_and_projection_have_independent_window_answers(seed,perspective):
    axis,angle = [1.,2.,-1.], 17.*seed-40.
    eye = [3.,-1.,2.]
    forward,up = rotate([0.,0.,-1.],axis,angle),rotate([0.,1.,0.],axis,angle)
    view = matrix.lookAt(V(eye),V([a+b for a,b in zip(eye,forward)]),V(up))
    camera = [.25,-.5,-(2.+seed)]
    world = [a+b for a,b in zip(rotate(camera,axis,angle),eye)]
    viewport = [11.,23.,320.,90.]
    projection = matrix.perspective(60.,2.,1.,20.) if perspective else matrix.orthographic(-2.,2.,-1.,1.,1.,20.)
    if perspective:
        cot = 1/math.tan(math.pi/6)
        ndc = [camera[0]*cot/(2*-camera[2]),camera[1]*cot/-camera[2],
               (21./19.*-camera[2]-40./19.)/-camera[2]]
    else:
        ndc = [camera[0]/2.,camera[1],(-2*camera[2]-21.)/19.]
    expected = [11.+(ndc[0]+1)*160.,23.+(ndc[1]+1)*45.,(ndc[2]+1)/2.]
    for model in (view,view.matrix):
        for proj in (projection,projection.matrix):
            assert matrix.project(V(world+[1.]),model,proj,viewport).vector == pytest.approx(expected,abs=2e-12)
            # Independently derived windows, not the implementation's projection result.
            assert matrix.unproject(*expected,model,proj,viewport).vector == pytest.approx(world,abs=2e-12)


@pytest.mark.parametrize('dimension', [2,3,4])
@pytest.mark.parametrize('t', [-.5,0.,.125,.5,1.,1.5])
def test_bezier_affine_vector_data_and_linear_interpolation(dimension,t):
    points = [[float((i+1)*(j+1)-3) for j in range(dimension)] for i in range(4)]
    rows = [[float(i == j)*2+(.5 if j == i+1 else 0.) for j in range(dimension+1)] for i in range(dimension+1)]
    rows[-1] = [float(j-2) for j in range(dimension)]+[1.]
    for degree,evaluator in ((2,bezier.quadraticBezierPoint),(3,bezier.cubicBezierPoint)):
        controls = [V(p) for p in points[:degree+1]]
        storages = [p.vector for p in controls]
        out = evaluator(t,*controls)
        expected = [float(x) for x in bernstein(points[:degree+1],t)]
        assert out.vector == pytest.approx(expected,abs=2e-14)
        transformed = [Vector(dimension).transform(p.vector,rows) for p in controls]
        after = evaluator(t,*transformed)
        assert after.vector == pytest.approx([float(x) for x in row(expected+[1.],rows)[:dimension]],abs=3e-14)
        assert all(p.vector is storage and p.vector == original for p,storage,original in zip(controls,storages,points))
        assert all(out.vector is not storage for storage in storages)
    interpolation = V(points[0])+((V(points[1])-V(points[0]))*t)
    assert interpolation.vector == [common.scalarLerp(a,b,t) for a,b in zip(points[0],points[1])]


@pytest.mark.parametrize('size', [2,3])
def test_adaptive_path_transformed_samples_and_repeated_ownership(size):
    controls = [V(p[:size]) for p in [[0.,0.,0.],[1.,2.,1.],[2.,-1.,2.],[3.,0.,3.]]]
    rows = matrix.Matrix(size+1).translate(V([2.,-1.] if size==2 else [2.,-1.,3.])).matrix
    moved = [Vector(size).transform(p.vector,rows) for p in controls]
    path = bezier.BezierPath(); path.setControlPoints(moved); path.minimum_sqr_distance = 1e-4
    samples = path.getDrawingPoints()[0]
    parameters = [(p.vector[0]-2.)/3. for p in samples]
    assert parameters == sorted(parameters) and parameters[0] == 0. and parameters[-1] == 1.
    for sample,t in zip(samples,parameters):
        expected = row([float(x) for x in bernstein([p.vector for p in controls],t)]+[1.],rows)
        assert sample.vector == pytest.approx([float(x) for x in expected[:size]],abs=3e-14)
    repeated = path.getDrawingPoints()[0]
    assert [p.vector for p in repeated] == [p.vector for p in samples]
    assert all(a is not b and a.vector is not b.vector for a,b in zip(repeated,samples))
    assert path.getControlPoints() is moved and all(s.vector is not c.vector for s in samples for c in moved)


@pytest.mark.parametrize('exponent', [-996,-500,0,500,996])
@pytest.mark.parametrize('size', [3,4])
def test_scaled_inverse_and_vector_application_remain_finite(exponent,size):
    base = [[3.,1.,0.],[0.,2.,1.],[1.,0.,2.]]
    if size == 4:
        base = [r+[0.] for r in base]+[[1.,-.5,.25,1.]]
    rows = [[math.ldexp(x,exponent) for x in r] for r in base]
    wrapped = matrix.Matrix(size,rows)
    inverted = wrapped.inverse()
    reference = inverse_reference(rows)
    for actual,expected in zip(inverted.matrix,reference):
        assert actual == pytest.approx([float(x) for x in expected],rel=3e-14,abs=0.)
    point = V([math.ldexp(float(i+1),-exponent) for i in range(size)])
    mapped = wrapped*point
    assert mapped.vector == pytest.approx([float(x) for x in row(point.vector,rows)],rel=3e-14,abs=3e-14)
    restored = inverted*mapped
    assert restored.vector == pytest.approx(point.vector,rel=4e-14,abs=0.)
    saved_export = wrapped.c_matrix
    assert wrapped.i_inverse() is wrapped and wrapped.c_matrix is not saved_export
    for actual,expected in zip(wrapped.c_matrix,inverted.c_matrix):
        assert list(actual) == list(expected)


@pytest.mark.parametrize('t', [0.,1e-12,.25,.5,.75,math.nextafter(1.,0.),1.])
@pytest.mark.parametrize('negate', [False,True])
def test_interpolated_rotations_feed_matrices_and_sh(t,negate):
    axis = [1.,2.,2.]
    q0 = quaternion.quat_from_axis_angle(axis,20.)
    q1 = quaternion.quat_from_axis_angle(axis,100.)
    if negate:
        q1 = q1.negate()
    coefficients = [[(i-4)*.1,(3-i)*.2,.05*i] for i in range(9)]
    for out,angle in ((quaternion.quat_slerp(q0,q1,t),20.+80.*t),
        (quaternion.squad4(q0,q1,quaternion.quat_from_axis_angle(axis,40.),
                           quaternion.quat_from_axis_angle(axis,80.),t),
         (1.-2*t*(1-t))*(20.+80.*t)+2*t*(1-t)*(40.+40.*t))):
        expected = rotate([.6,0.,.8],axis,angle)
        assert (out.toMatrix()*V([.6,0.,.8,0.])).vector == pytest.approx(expected+[0.],abs=3e-14)
        rotated = sh.rotate_coefficients(coefficients,out)
        direction = [.3,.4,math.sqrt(.75)]
        assert sh.reconstruct(rotated,direction) == pytest.approx(evaluate(coefficients,rotate(direction,axis,-angle)),abs=3e-14)
        convolved = sh.convolve_diffuse(rotated)
        factors = [math.pi]+[2*math.pi/3]*3+[math.pi/4]*5
        reference_coefficients = [[x*f for x in c] for c,f in zip(coefficients,factors)]
        assert sh.reconstruct(convolved,direction) == pytest.approx(
            evaluate(reference_coefficients,rotate(direction,axis,-angle)),abs=4e-14)
        assert all(r is not c for r in convolved for c in rotated)


@pytest.mark.parametrize('component', [1e-9,1e-100,1e-300,math.ulp(0.)])
@pytest.mark.parametrize('pole', [-1.,1.])
def test_near_pole_matrix_quaternion_and_sh_consistency(component,pole):
    direction = [component,0.,pole*math.sqrt(1-component*component)]
    q = quaternion.Quaternion([0.,1.,0.,0.])  # Exact half-turn: (x,y,z)->(x,-y,-z).
    expected = [direction[0],-direction[1],-direction[2]]
    assert quaternion.quat_rotate_vector(q,V(direction)).vector == expected
    assert (q.toMatrix()*V(direction+[0.])).vector == expected+[0.]
    coefficients = [[0.]*3 for _ in range(9)]; coefficients[3] = [1.,2.,-1.]
    out = sh.rotate_coefficients(coefficients,q)
    actual = sh.reconstruct(out,V(expected))
    reference = evaluate(coefficients,direction)
    for a,b in zip(actual,reference):
        assert abs(a-b) <= 8*math.ulp(b)


@pytest.mark.parametrize('x', [-.8,-.25,0.,.25,.8])
def test_legendre_phase_and_sh_normalization_against_polynomials(x):
    # Explicit low-order Rodrigues derivatives; higher orders use existing independent audits.
    theta,phi = math.acos(x),.37
    polynomials = [1.,x,(3*x*x-1)/2.,(5*x*x*x-3*x)/2.]
    for l,expected in enumerate(polynomials):
        p = legendre.Legendre(l,0,x)
        assert p.run() == pytest.approx(expected,abs=3e-15)
        assert sh.SPH(l,0,theta,phi) == pytest.approx(math.sqrt((2*l+1)/(4*math.pi))*expected,abs=3e-15)
        assert (p.P,p.PM1,p.PML) == (1.,0.,0.)
    assert legendre.Legendre(2,1,x).run() == pytest.approx(-3*x*math.sqrt(1-x*x),abs=3e-15)
    direction = [math.sqrt(1-x*x)*math.cos(phi),math.sqrt(1-x*x)*math.sin(phi),x]
    reference = basis(direction)
    assert [sh.SPH(l,m,theta,phi) for l in range(3) for m in range(-l,l+1)] == pytest.approx(reference,abs=3e-15)


@pytest.mark.parametrize('seed',range(5))
def test_independent_axis_cubature_projects_rotates_and_convolves_rgb(seed):
    directions = [[s if i==axis else 0. for i in range(3)] for axis in range(3) for s in (-1.,1.)]
    samples = []
    colors = []
    for d in directions:
        sample = sh.SPHSample(math.acos(d[2]),math.atan2(d[1],d[0]),V(d),9)
        sample.values = basis(d)  # Independent Cartesian values, not gem SPH.
        samples.append(sample)
        colors.append([2.+.5*d[0],3.-.2*d[1],4.+.75*d[2]])
    coefficients = sh.project_radiance(samples,colors)
    expected = [[0.]*3 for _ in range(9)]
    expected[0] = [x*math.sqrt(4*math.pi) for x in (2.,3.,4.)]
    k = math.sqrt(4*math.pi/3)
    expected[1][1],expected[2][2],expected[3][0] = .2*k,.75*k,-.5*k
    for actual,reference in zip(coefficients,expected):
        assert actual == pytest.approx(reference,abs=5e-15)
    assert sh.project_radiance(samples,colors,[4*math.pi/6]*6) == coefficients
    axis,angle = [1.,-2.,2.],23.*seed-37.
    q = quaternion.quat_from_axis_angle(axis,angle)
    irradiance = sh.convolve_diffuse(sh.rotate_coefficients(coefficients,q))
    d = [.3,.4,math.sqrt(.75)]
    inverse_d = rotate(d,axis,-angle)
    expected_value = [2*math.pi+(2*math.pi/3)*.5*inverse_d[0],
                      3*math.pi-(2*math.pi/3)*.2*inverse_d[1],
                      4*math.pi+(2*math.pi/3)*.75*inverse_d[2]]
    assert sh.reconstruct(irradiance,V(d)) == pytest.approx(expected_value,abs=2e-14)
    assert colors == [[2.+.5*d[0],3.-.2*d[1],4.+.75*d[2]] for d in directions]
    assert all(c is not e for c in coefficients for e in colors)


@pytest.mark.parametrize('width,height',[(3,5),(6,2)])
def test_legacy_raw_probe_requires_explicit_conversion_before_rotation(tmp_path,width,height):
    import struct
    pixels = [[[1.+i/8.,.5+j/4.,.25] for i in range(width)] for j in range(height)]
    raw = tmp_path/'probe.float32'
    raw.write_bytes(struct.pack('%df'%(3*width*height),*(v for row in pixels for p in row for v in p)))
    legacy = sh.SPH_IrradianceMapCoeff(str(raw),width,height)
    canonical = sh.legacy_to_canonical(legacy.coeffs)
    independent = [[0.]*3 for _ in range(9)]
    for j,row_pixels in enumerate(pixels):
        for i,color in enumerate(row_pixels):
            u,v = 2*(i+.5)/width-1,1-2*(j+.5)/height
            r = math.hypot(u,v)
            if r > 1: continue
            theta = math.pi*r
            direction = [math.sin(theta)*u/r,math.sin(theta)*v/r,math.cos(theta)] if r else [0.,0.,1.]
            weight = 4*math.pi**2/(width*height)*(math.sin(theta)/theta if theta else 1.)
            for k,b in enumerate(basis(direction)):
                for channel in range(3): independent[k][channel] += color[channel]*weight*b
    for actual,expected in zip(canonical,independent):
        assert actual == pytest.approx(expected,abs=3e-14)
    for actual,expected in zip(sh.project_angular_probe(pixels),independent):
        assert actual == pytest.approx(expected,abs=3e-14)
    q = quaternion.quat_from_axis_angle([0.,0.,1.],90.)
    result = sh.rotate_coefficients(canonical,q)
    assert sh.reconstruct(result,[0.,1.,0.]) == pytest.approx(evaluate(independent,[1.,0.,0.]),abs=4e-14)
    old = [c[:] for c in legacy.coeffs]
    assert legacy.calculateCoefficients() is None and legacy.coeffs == old
    assert legacy.hdr == pixels


def test_shared_storage_across_vector_quaternion_and_ctypes_boundaries():
    data = [1.,2.,2.,0.]
    v,q = Vector(4,data),quaternion.Quaternion(data)
    assert q.data is v.vector is data
    result = v.normalize()
    assert result.vector == [1/3.,2/3.,2/3.,0.] and q.data is data
    assert q.i_conjugate() is q and q.data is not data
    assert data == v.vector == [1.,2.,2.,0.]
    assert v.i_clamp(4,data,[0.]*4,[1.]*4) is v and v.vector is not data
    assert data == [1.,2.,2.,0.] and q.data == [1.,-2.,-2.,-0.]
    snapshot = common.conv_list(data,common.GLfloat)
    data[0] = 9.
    assert snapshot[0] == 1. and result.vector[0] == 1/3.


def test_direct_zero_fallbacks_do_not_hide_geometric_degeneracy():
    assert Vector(3).normalize().vector == [0.]*3
    assert quaternion.Quaternion([0.]*4).normalize().data == [1.,0.,0.,0.]
    cases = [lambda: ray.Ray(V([0.]*3),Vector(3)),
             lambda: plane.Plane().normalize(),
             lambda: plane.Plane().fromPoints(V([0.]*3),V([1.]*3),V([2.]*3)),
             lambda: quaternion.quat_from_axis_angle(Vector(3),45.),
             lambda: matrix.Matrix(4).rotate(Vector(3),45.),
             lambda: matrix.lookAt(Vector(3),Vector(3),V([0.,1.,0.]))]
    for operation in cases:
        with pytest.raises(ZeroDivisionError): operation()
    with pytest.raises(ValueError): sh.rotate_coefficients([1.],quaternion.Quaternion([0.]*4))
    with pytest.raises(ValueError): sh.reconstruct([[1.]*3],Vector(3))


@pytest.mark.parametrize('api', ['free_slerp','method_slerp','squad4'])
@pytest.mark.parametrize('sign', [-1.,1.])
@pytest.mark.parametrize('representation', [-1.,1.])
@pytest.mark.defect('4G5-A01: subnormal SLERP separation loses the endpoint rotation')
def test_subnormal_interpolation_endpoint_retains_representable_rotation(api,sign,representation):
    tiny = math.ulp(0.)
    constructed = quaternion.quat_from_axis_angle([1.,0.,0.],math.degrees(2*sign*tiny))
    assert constructed.data[1] == sign*tiny and math.hypot(*constructed.data) == 1.
    endpoint = constructed if representation == 1. else constructed.negate()
    start = quaternion.Quaternion()
    storage,saved = endpoint.data,endpoint.data[:]
    if api == 'free_slerp':
        out = quaternion.quat_slerp(start,endpoint,1.)
    elif api == 'method_slerp':
        out = start.slerp(endpoint,1.)
    else:
        out = quaternion.squad4(start,endpoint,start,endpoint,1.)
    assert endpoint.data is storage and endpoint.data == saved
    assert out is not endpoint and out.data is not storage
    # At t=1, q*[0,Y]*q* has Z=2*w*x exactly; Fraction supplies the answer.
    expected = [0.,1.,float(2*F(endpoint.data[0])*F(endpoint.data[1]))]
    assert quaternion.quat_rotate_vector(endpoint,V([0.,1.,0.])).vector == expected
    assert quaternion.quat_rotate_vector(out,V([0.,1.,0.])).vector == expected
