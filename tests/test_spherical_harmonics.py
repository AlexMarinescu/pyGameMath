import math
import random
import struct
import pytest
from gem import spherical_harmonics as sh
from gem.vector import Vector

def test_generate_samples(monkeypatch):
    monkeypatch.setattr(random,'random',lambda:0.5)
    samples = sh.GenerateSamples(2,3)
    assert len(samples) == 4
    for sample in samples:
        assert sample.dir.magnitude() == pytest.approx(1)
        assert len(sample.values) == 9
        assert sample.values[0] == pytest.approx(1/math.sqrt(4*math.pi))


def test_rectangular_irradiance_file(tmp_path):
    path = tmp_path/'probe.float'
    path.write_bytes(struct.pack('12f',*([1.0]*12)))
    result = sh.SPH_IrradianceMapCoeff(str(path),2,1)
    assert len(result.hdr) == 1
    assert len(result.hdr[0]) == 2


def polynomial_basis(x,y,z):
    # Independent Cartesian formulas with exact normalization constants.
    a=math.sqrt(3/(4*math.pi));b=math.sqrt(15/(4*math.pi))
    return [1/math.sqrt(4*math.pi),-a*y,a*z,-a*x,b*x*y,-b*y*z,
            math.sqrt(5/(16*math.pi))*(3*z*z-1),-b*x*z,
            math.sqrt(15/(16*math.pi))*(x*x-y*y)]


def sphere_grid(nz=40,nphi=80):
    samples=[];directions=[]
    for i in range(nz):
        z=-1+2*(i+.5)/nz
        radius=math.sqrt(1-z*z)
        for j in range(nphi):
            phi=2*math.pi*(j+.5)/nphi
            direction=[radius*math.cos(phi),radius*math.sin(phi),z]
            sample=sh.SPHSample(math.acos(z),phi,Vector(3,direction),9)
            sample.values=polynomial_basis(*direction)
            samples.append(sample);directions.append(direction)
    return samples,directions


@pytest.mark.parametrize('direction',[(1,0,0),(0,1,0),(0,0,1),(0,0,-1),(.3,.4,math.sqrt(.75))])
def test_basis_index_sign_and_cartesian_reference(direction):
    x,y,z=direction
    actual=[sh.SPH(l,m,math.acos(z),math.atan2(y,x)) for l in range(3) for m in range(-l,l+1)]
    assert actual==pytest.approx(polynomial_basis(x,y,z),abs=2e-15)
    assert [l*(l+1)+m for l in range(3) for m in range(-l,l+1)]==list(range(9))


def test_constant_radiance_coefficients_and_constant_diffuse():
    samples,directions=sphere_grid()
    color=[2.,3.,5.]
    coeff=sh.project_radiance(samples,[color]*len(samples))
    assert coeff[0]==pytest.approx([v*math.sqrt(4*math.pi) for v in color])
    # Midpoint integration of the degree-two z polynomial has O(nz^-2) error.
    assert coeff[6]==pytest.approx([-math.sqrt(5*math.pi)*v/40**2 for v in color],abs=1e-14)
    assert max(abs(v) for i,row in enumerate(coeff) if i not in (0,6) for v in row)<1e-14
    exact=[[0.]*3 for _ in range(9)];exact[0]=[v*math.sqrt(4*math.pi) for v in color]
    irradiance=sh.convolve_diffuse(exact)
    for direction in [(1,0,0),(0,1,0),(0,0,1)]:
        assert sh.reconstruct(exact,direction)==pytest.approx(color)
        assert sh.reconstruct(irradiance,direction)==pytest.approx([math.pi*v for v in color])
    assert exact[0]==[v*math.sqrt(4*math.pi) for v in color]
    assert all(a is not b for a,b in zip(exact,irradiance))


def test_asymmetric_rgb_projection_orientation_and_reconstruction():
    samples,directions=sphere_grid(80,80)
    # Analytic integrals: int x^2 dOmega=int y^2=int z^2=4pi/3.
    radiance=[[1+2*x,2+3*y,3+4*z] for x,y,z in directions]
    coeff=sh.project_radiance(samples,radiance)
    expected=[[0.]*3 for _ in range(9)]
    expected[0]=[math.sqrt(4*math.pi)*v for v in [1,2,3]]
    expected[3][0]=-2*math.sqrt(4*math.pi/3)
    expected[1][1]=-3*math.sqrt(4*math.pi/3)
    expected[2][2]=4*math.sqrt(4*math.pi/3)
    for actual,reference in zip(coeff,expected):
        assert actual==pytest.approx(reference,rel=4e-4,abs=.003)
    for direction in [(1,0,0),(0,1,0),(0,0,1),(-1,0,0)]:
        x,y,z=direction
        assert sh.reconstruct(expected,direction)==pytest.approx([1+2*x,2+3*y,3+4*z])
        assert sh.reconstruct(sh.convolve_diffuse(expected),direction)==pytest.approx(
            [math.pi+4*math.pi*x/3,2*math.pi+2*math.pi*y,3*math.pi+8*math.pi*z/3])


def test_explicit_weights_ownership_and_repeated_projection():
    samples,directions=sphere_grid(2,4)
    radiance=[[i+1,2*i+1,3*i+1] for i in range(len(samples))]
    weights=[.1*(i+1) for i in range(len(samples))]
    before=[s.values[:] for s in samples];vectors=[s.dir.vector[:] for s in samples]
    colors=[c[:] for c in radiance];savedweights=weights[:]
    expected=[[sum(c[channel]*w*polynomial_basis(*d)[index] for c,w,d in zip(radiance,weights,directions))
               for channel in range(3)] for index in range(9)]
    for _ in range(2):
        actual=sh.project_radiance(samples,radiance,weights)
        for a,b in zip(actual,expected):assert a==pytest.approx(b)
    assert [s.values for s in samples]==before
    assert [s.dir.vector for s in samples]==vectors
    assert radiance==colors and weights==savedweights


def test_sampling_seeded_repeatability_storage_and_solid_angle_strata(capsys):
    random.seed(127);a=sh.GenerateSamples(8,3)
    random.seed(127);b=sh.GenerateSamples(8,3)
    assert [(p.theta,p.phi,p.dir.vector,p.values) for p in a]==[(p.theta,p.phi,p.dir.vector,p.values) for p in b]
    assert capsys.readouterr().out==''
    for index,sample in enumerate(a):
        i,j=divmod(index,8)
        u=(1-sample.dir.vector[2])/2;v=sample.phi/(2*math.pi)
        assert i/8<=u<(i+1)/8 and j/8<=v<(j+1)/8
        assert sample.dir.magnitude()==pytest.approx(1,abs=2e-15)
        assert sample.values==pytest.approx(polynomial_basis(*sample.dir.vector),abs=2e-15)
        assert not hasattr(sample.dir,'vec')
    assert len({id(s.dir.vector) for s in a})==64
    supplied=Vector(3,[1,0,0]);sample=sh.SPHSample(0,0,supplied,9)
    assert sample.dir is supplied # historical constructor ownership


def angular_image(width,height,function):
    image=[]
    for row in range(height):
        values=[]
        v=1-2*(row+.5)/height
        for col in range(width):
            u=2*(col+.5)/width-1;r=math.sqrt(u*u+v*v)
            theta=math.pi*r;phi=math.atan2(v,u)
            d=(math.sin(theta)*math.cos(phi),math.sin(theta)*math.sin(phi),math.cos(theta))
            values.append(list(function(*d)))
        image.append(values)
    return image


@pytest.mark.parametrize('width,height',[(80,80),(120,60)])
def test_angular_probe_constant_weight_asymmetry_and_rgb(width,height):
    image=angular_image(width,height,lambda x,y,z:[1+x,2+y,3+z])
    before=[[pixel[:] for pixel in row] for row in image]
    coeff=sh.project_angular_probe(image)
    expected=[[0.]*3 for _ in range(9)]
    expected[0]=[math.sqrt(4*math.pi)*i for i in [1,2,3]]
    expected[3][0]=expected[1][1]=-math.sqrt(4*math.pi/3)
    expected[2][2]=math.sqrt(4*math.pi/3)
    for a,b in zip(coeff,expected):assert a==pytest.approx(b,abs=.004)
    assert image==before
    assert sh.project_angular_probe(image)==coeff


@pytest.mark.parametrize('row,col',[(1,2),(0,1),(1,1),(2,1),(1,0)])
def test_pixel_centers_orientation_and_solid_angle(row,col):
    image=[[[0.]*3 for _ in range(3)] for _ in range(3)]
    image[row][col]=[1,2,3]
    u=2*(col+.5)/3-1;v=1-2*(row+.5)/3;r=math.hypot(u,v)
    theta=math.pi*r
    direction=(math.sin(theta)*u/r,math.sin(theta)*v/r,math.cos(theta)) if r else (0,0,1)
    weight=4*math.pi**2/9*(math.sin(theta)/theta if theta else 1)
    actual=sh.project_angular_probe(image)
    for a,value in zip(actual,polynomial_basis(*direction)):
        assert a==pytest.approx([value*weight*i for i in [1,2,3]],abs=2e-14)
    if col==2:assert actual[3][0]<0
    if row==0:assert actual[1][0]<0


def test_outside_disk_ignored_and_no_weight_renormalization():
    image=[[[0.]*3 for _ in range(4)] for _ in range(4)];image[0][0]=[1,1,1]
    assert sh.project_angular_probe(image)==[[0.]*3 for _ in range(9)]
    # One center pixel has midpoint quadrature area 4pi^2, not 4pi.
    coeff=sh.project_angular_probe([[[1,1,1]]],1)
    assert coeff[0]==pytest.approx([4*math.pi**2/math.sqrt(4*math.pi)]*3)


def test_legacy_loader_orientation_conversion_and_recalculation(tmp_path):
    image=angular_image(12,8,lambda x,y,z:[1+x,2+y,3+z])
    values=[channel for row in image for pixel in row for channel in pixel]
    file=tmp_path/'probe.float';file.write_bytes(struct.pack('%df'%len(values),*values))
    p=sh.SPH_IrradianceMapCoeff(str(file),12,8)
    assert len(p.hdr)==8 and len(p.hdr[0])==12
    for a,b in zip([c for row in p.hdr for pixel in row for c in pixel],values):assert a==pytest.approx(b,rel=1e-6)
    original=[c[:] for c in p.coeffs]
    p.calculateCoefficients();assert p.coeffs==original
    canonical=sh.project_angular_probe(p.hdr)
    converted=sh.legacy_to_canonical(p.coeffs)
    for a,b in zip(converted,canonical):assert a==pytest.approx(b,abs=1e-13)
    assert p.coeffs==original
    for a,b in zip(converted,p.coeffs):assert a is not b
    p.updateCoefficients([1,2,3],1,1,0,0)
    assert p.coeffs[3][0]==pytest.approx(original[3][0]+.488603)


def test_diffuse_factors_no_double_convolution_and_no_mutation():
    values=[[1.,2.,3.] for _ in range(9)];before=[v[:] for v in values]
    irradiance=sh.convolve_diffuse(values)
    for i,row in enumerate(irradiance):
        factor=[math.pi,2*math.pi/3,math.pi/4][math.isqrt(i)]
        assert row==pytest.approx([factor,2*factor,3*factor])
    direction=Vector(3,[0,0,1]);storage=direction.vector
    actual=sh.reconstruct(irradiance,direction)
    expected=[sum(row[c]*basis for row,basis in zip(irradiance,polynomial_basis(0,0,1))) for c in range(3)]
    assert actual==pytest.approx(expected)
    assert direction.vector is storage and storage==[0,0,1]
    assert values==before


def test_orthogonality_normalization_and_midpoint_weight_sum():
    samples,_=sphere_grid(120,48)
    weight=4*math.pi/len(samples)
    assert weight*len(samples)==pytest.approx(4*math.pi)
    for i in range(9):
        for j in range(i+1):
            integral=math.fsum(s.values[i]*s.values[j]*weight for s in samples)
            assert integral==pytest.approx(float(i==j),abs=.0004)


@pytest.mark.parametrize('call',[
    lambda:sh.project_radiance([],[]),
    lambda:sh.project_radiance([sh.SPHSample(0,0,Vector(3,[0,0,1]),9)],[]),
    lambda:sh.project_radiance([sh.SPHSample(0,0,Vector(3,[0,0,1]),2)],[[1,1,1]]),
    lambda:sh.project_radiance([sh.SPHSample(0,0,Vector(3,[0,0,1]),1)],[[1,1,1]],[-1]),
    lambda:sh.project_angular_probe([]),
    lambda:sh.project_angular_probe([[[1,1,1]],[]]),
    lambda:sh.project_angular_probe([[[1,float('nan'),1]]]),
    lambda:sh.project_angular_probe([[[1,1,1]]],0),
    lambda:sh.reconstruct([[1,1,1]],[0,0,0]),
    lambda:sh.reconstruct([[1,1,1]],[0,0,2]),
    lambda:sh.reconstruct([[1,1,1]],[0,float('inf'),1]),
    lambda:sh.reconstruct([[1,1,1]]*2,[0,0,1]),
    lambda:sh.convolve_diffuse([[1,1,1]]*16),
    lambda:sh.legacy_to_canonical([[1,1,1]]*4),
])
def test_new_api_invalid_inputs(call):
    with pytest.raises(ValueError):call()


def test_black_environment_zero_weights_and_signed_hdr_inputs():
    samples,_=sphere_grid(2,4)
    assert sh.project_radiance(samples,[[1,1,1]]*8,[0]*8)==[[0.]*3 for _ in range(9)]
    assert sh.project_angular_probe([[[0,0,0]]])==[[0.]*3 for _ in range(9)]
    # No clipping policy is introduced for mathematical RGB input values.
    assert sh.project_radiance(samples,[[-1,2,0]]*8)[0]==pytest.approx([-math.sqrt(4*math.pi),2*math.sqrt(4*math.pi),0])


@pytest.mark.parametrize('width,height',[(0,1),(-1,1),(1,0),(1.5,2)])
def test_invalid_raw_dimensions(tmp_path,width,height):
    f=tmp_path/'short.float';f.write_bytes(b'')
    with pytest.raises(ValueError):sh.SPH_IrradianceMapCoeff(str(f),width,height)


def test_truncated_raw_probe(tmp_path):
    f=tmp_path/'short.float';f.write_bytes(struct.pack('f',1))
    with pytest.raises(ValueError,match='enough RGB'):sh.SPH_IrradianceMapCoeff(str(f),2,1)


def test_compatibility_import_identity_and_input_types():
    from gem.experimental import sph,sph_sample,sph_irradiance_map
    assert sph.Factorial is sh.Factorial and sph.K is sh.K and sph.SPH is sh.SPH
    assert sph.Legendre is sh.Legendre
    assert sph_sample.SPHSample is sh.SPHSample
    assert sph_sample.GenerateSamples is sh.GenerateSamples
    assert sph_irradiance_map.SPH_IrradianceMapCoeff is sh.SPH_IrradianceMapCoeff


def test_constant_hdr_to_diffuse_convergence():
    errors=[]
    for width in [16,32,64,128]:
        image=[[[1.,2.,3.] for _ in range(width)] for _ in range(width//2)]
        radiance=sh.project_angular_probe(image)
        irradiance=sh.convolve_diffuse(radiance)
        errors.append(max(abs(value-math.pi*(channel+1))
                          for direction in [(1,0,0),(0,1,0),(0,0,1)]
                          for channel,value in enumerate(sh.reconstruct(irradiance,direction))))
    assert errors==sorted(errors,reverse=True)
    assert errors[-1]<5e-5


def test_finite_mixed_scale_rgb_and_sample_storage():
    samples,_=sphere_grid(2,4)
    storage=[(sample.values,sample.dir,sample.dir.vector) for sample in samples]
    colors=[[1e200,1e-200,1.]]*len(samples)
    coefficients=sh.project_radiance(samples,colors)
    assert coefficients[0]==pytest.approx([math.sqrt(4*math.pi)*v for v in colors[0]],rel=2e-15,abs=0)
    for sample,(values,direction,components) in zip(samples,storage):
        assert sample.values is values and sample.dir is direction and sample.dir.vector is components
