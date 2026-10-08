import json
import math
import re
import struct
import zlib
from pathlib import Path
import pytest
from gem import spherical_harmonics as sh
from gem.quaternion import Quaternion
from examples.hdr_sh import environment as env, reference, visualize

ROOT=Path(__file__).resolve().parents[1]


def shader_evaluate(coefficients,direction):
    # Evaluate the actual checked-in GLSL expression, not a second hand-written
    # SH implementation. Each operation is rounded to binary32 on the CPU.
    class Float32(float):
        def __new__(cls,v):return float.__new__(cls,struct.unpack('f',struct.pack('f',v))[0])
        def __add__(self,v):return Float32(float(self)+float(v))
        __radd__=__add__
        def __sub__(self,v):return Float32(float(self)-float(v))
        def __rsub__(self,v):return Float32(float(v)-float(self))
        def __mul__(self,v):return Float32(float(self)*float(v))
        __rmul__=__mul__
    source=(ROOT/'examples/hdr_sh/diffuse.glsl').read_text()
    expression=re.search(r'return (irradianceSH\[0\].*?);',source,re.S).group(1)
    result=[]
    for channel in range(3):
        expr=re.sub(r'irradianceSH\[(\d+)\]',lambda m:'c['+m.group(1)+']',expression)
        x,y,z=map(Float32,direction);c=[Float32(row[channel]) for row in coefficients]
        result.append(float(eval('('+expr+')',{'__builtins__':{}},{'c':c,'x':x,'y':y,'z':z})))
    return result


def test_latlong_orientation_and_exact_pixel_solid_angle():
    theta,phi,direction,weight=env.latlong_pixel(0,0,4,2)
    assert theta==math.pi/4 and phi==math.pi/4
    assert direction==pytest.approx([.5,.5,math.sqrt(.5)])
    assert weight==pytest.approx(math.pi/2)
    assert env.latlong_pixel(0,1,4,2)[2][0]<0
    assert env.latlong_pixel(1,0,4,2)[2][2]<0
    assert sum(env.latlong_pixel(r,c,13,7)[3] for r in range(7) for c in range(13))==pytest.approx(4*math.pi)


def test_constant_linear_hdr_to_lambertian_and_display():
    image=env.synthetic_environment(64,32,[2,3,4],[[0,0,0]]*3)
    result=reference.run(image,'latlong',Quaternion(),[.2,.4,.6])
    assert result['radiance_coefficients'][0]==pytest.approx([v*math.sqrt(4*math.pi) for v in [2,3,4]])
    for entry in result['evaluations']:
        assert entry['irradiance_linear_rgb']==pytest.approx([math.pi*v for v in [2,3,4]],abs=.009)
        assert entry['lambertian_linear_rgb']==pytest.approx([.4,1.2,2.4],abs=.002)
    assert reference.linear_to_srgb(.003)==pytest.approx(.03876)
    assert reference.linear_to_srgb(.5)==pytest.approx(.7353569830524495)
    assert reference.linear_to_srgb(2)>1 # no premature clipping/tone mapping


def test_fixture_axis_rotation_analytical_reference_repeatability_and_ownership():
    image=env.load_fixture(ROOT/'examples/hdr_sh/fixtures/asymmetric.json')
    before=json.dumps(image)
    orientation=Quaternion([math.sqrt(.5),0,0,math.sqrt(.5)]);storage=orientation.data
    a=reference.run(image,'latlong',orientation,[.5]*3)
    b=reference.run(image,'latlong',orientation,[.5]*3)
    assert a==b and json.dumps(image)==before and orientation.data is storage
    for entry in a['evaluations']:
        x,y,z=entry['normal']
        # Rotated [2+X,3+Y,4+Z] becomes [2+Y,3-X,4+Z].
        expected=[math.pi*2+2*math.pi*y/3,math.pi*3-2*math.pi*x/3,math.pi*4+2*math.pi*z/3]
        assert entry['irradiance_linear_rgb']==pytest.approx(expected,abs=.013)
        inverse=[y,-x,z]
        original=sh.convolve_diffuse(a['radiance_coefficients'])
        assert entry['irradiance_linear_rgb']==pytest.approx(sh.reconstruct(original,inverse),abs=2e-14)
        assert shader_evaluate(a['irradiance_coefficients'],entry['normal'])==pytest.approx(entry['irradiance_linear_rgb'],rel=2e-6,abs=2e-6)
    exported=reference.glsl_uniforms(a['irradiance_coefficients'])
    assert exported.count('vec3(')==9


def test_single_pixel_direction_and_rgb_independence():
    image=[[[0.]*3 for _ in range(4)] for _ in range(2)];image[0][0]=[1,2,4]
    coefficients=env.project_latlong(image)
    x,y,z=.5,.5,math.sqrt(.5)
    a=math.sqrt(3/(4*math.pi));b=math.sqrt(15/(4*math.pi))
    basis=[1/math.sqrt(4*math.pi),-a*y,a*z,-a*x,b*x*y,-b*y*z,
           math.sqrt(5/(16*math.pi))*(3*z*z-1),-b*x*z,math.sqrt(15/(16*math.pi))*(x*x-y*y)]
    for row,value in zip(coefficients,basis):assert row==pytest.approx([math.pi/2*value*c for c in [1,2,4]],abs=2e-15)
    assert coefficients[1][0]<0 and coefficients[3][0]<0


@pytest.mark.parametrize('rle',[False,True])
def test_actual_rgbe_decoding(tmp_path,rle):
    width=8 if rle else 4
    header=b'#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y 2 +X '+str(width).encode()+b'\n'
    if rle:row=b'\x02\x02\x00\x08'+b''.join(bytes([136,v]) for v in [32,64,128,131])
    else:row=bytes([32,64,128,131])*width
    file=tmp_path/'fixture.hdr';file.write_bytes(header+row*2)
    assert env.load_rgbe(file)==[[[1.,2.,4.] for _ in range(width)] for _ in range(2)]


def test_raw_angular_adapter_matches_core_and_rejects_wrong_fixture_mapping(tmp_path):
    file=tmp_path/'probe.float';file.write_bytes(struct.pack('6f',1,2,4,1,2,4))
    image=env.load_raw(file,2,1)
    result=reference.run(image,'angular',Quaternion(),[.5]*3)
    assert result['radiance_coefficients']==sh.project_angular_probe(image)
    with pytest.raises(ValueError):reference.run(image,'wrong',Quaternion(),[.5]*3)


def decode_png(file):
    data=Path(file).read_bytes();assert data[:8]==b'\x89PNG\r\n\x1a\n'
    position=8;compressed=b'';width=height=None
    while position<len(data):
        size=struct.unpack('>I',data[position:position+4])[0];kind=data[position+4:position+8]
        payload=data[position+8:position+8+size]
        crc=struct.unpack('>I',data[position+8+size:position+12+size])[0]
        assert crc==zlib.crc32(kind+payload)&0xffffffff
        if kind==b'IHDR':width,height,depth,color,_,_,_=struct.unpack('>IIBBBBB',payload);assert (depth,color)==(8,2)
        if kind==b'IDAT':compressed+=payload
        position+=12+size
    raw=zlib.decompress(compressed);stride=width*3+1
    assert len(raw)==height*stride
    assert all(raw[r*stride]==0 for r in range(height))
    return width,height,b''.join(raw[r*stride+1:(r+1)*stride] for r in range(height))


def test_cpu_visualization_known_rotation_and_committed_outputs(tmp_path):
    image=visualize.directional_environment(64,32)
    assert max(v for row in image for pixel in row for v in pixel)>1
    c=env.project_latlong(image);q=Quaternion([math.sqrt(.5),0,0,math.sqrt(.5)])
    a,sa=visualize.render_sphere(sh.convolve_diffuse(c),48)
    b,sb=visualize.render_sphere(sh.convolve_diffuse(sh.rotate_coefficients(c,q)),48)
    assert a!=b and sa['maximum_byte']==sb['maximum_byte']
    assert sa['linear_luminance_centroid_pixels'][0]>24
    assert sb['linear_luminance_centroid_pixels'][1]<24
    file=tmp_path/'sphere.png';visualize.write_png(file,48,48,a)
    assert decode_png(file)==(48,48,a)
    originals=[decode_png(ROOT/'examples/output'/name) for name in ['sh_original.png','sh_rotated.png']]
    assert originals[0][:2]==originals[1][:2]==(192,192)
    assert originals[0][2]!=originals[1][2]
    assert sum(originals[0][2])==sum(originals[1][2])


def test_glsl_lambertian_factor_is_applied_once():
    source=(ROOT/'examples/hdr_sh/diffuse.glsl').read_text()
    factor=float(re.search(r'\* ([0-9.]+); // 1/pi',source).group(1))
    assert factor==pytest.approx(1/math.pi,abs=1e-16)
    assert reference.reflected_radiance([math.pi,2*math.pi,3*math.pi],[.2,.4,.6])==pytest.approx([.2,.8,1.8])


@pytest.mark.parametrize('payload',[
    b'not HDR',
    b'#?RADIANCE\nFORMAT=32-bit_rle_xyze\n\n-Y 1 +X 1\n'+bytes(4),
    b'#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n+Y 1 +X 1\n'+bytes(4),
    b'#?RADIANCE\nFORMAT=32-bit_rle_rgbe\n\n-Y 1 +X 8\n\x02\x02\x00\x08\x00',
])
def test_unsupported_and_malformed_rgbe(tmp_path,payload):
    file=tmp_path/'invalid.hdr';file.write_bytes(payload)
    with pytest.raises(ValueError):env.load_rgbe(file)


@pytest.mark.parametrize('kind',['L0','L1','rotated-L1','red-only'])
def test_controlled_render_values_before_tone_mapping_and_pixel_bytes(kind):
    norm=math.sqrt(4*math.pi);a=math.sqrt(3/(4*math.pi))
    coefficients=[[norm*v for v in [1.,2.,3.]]]+[[0.,0.,0.] for _ in range(3)]
    if kind in ('L1','rotated-L1'):coefficients[3][0]=-.4/a
    if kind=='rotated-L1':coefficients=sh.rotate_coefficients(coefficients,Quaternion([math.sqrt(.5),0,0,math.sqrt(.5)]))
    if kind=='red-only':coefficients=[[norm,0,0]]
    pixels,stats=visualize.render_sphere(coefficients,48,albedo=(.5,.5,.5),exposure=.7)
    for entry in stats['reference_pixels']:
        x,y,z=entry['normal']
        expected=[1.,2.,3.]
        if kind=='L1':expected[0]+=.4*x
        if kind=='rotated-L1':expected[0]+=.4*y
        if kind=='red-only':expected=[1.,0.,0.]
        assert entry['irradiance_linear_rgb']==pytest.approx(expected,abs=1e-14)
        reflected=[v*.5/math.pi for v in expected]
        assert entry['lambertian_linear_rgb']==pytest.approx(reflected,abs=1e-14)
        display=[]
        for value in reflected:
            mapped=value*.7/(1+value*.7)
            srgb=12.92*mapped if mapped<=.0031308 else 1.055*mapped**(1/2.4)-.055
            display.append(round(srgb*255))
        col,row=entry['pixel'];offset=(row*48+col)*3
        assert list(pixels[offset:offset+3])==display


def test_render_changes_with_coefficients_without_channel_bleed():
    base=[[math.sqrt(4*math.pi),math.sqrt(4*math.pi),math.sqrt(4*math.pi)]]
    changed=[[2*base[0][0],base[0][1],base[0][2]]]
    a,sa=visualize.render_sphere(base,48)
    b,sb=visualize.render_sphere(changed,48)
    assert a!=b
    assert a[1::3]==b[1::3] and a[2::3]==b[2::3]
    for x,y in zip(sa['reference_pixels'],sb['reference_pixels']):
        assert y['irradiance_linear_rgb'][0]==pytest.approx(2*x['irradiance_linear_rgb'][0])
        assert y['irradiance_linear_rgb'][1:]==x['irradiance_linear_rgb'][1:]


def test_committed_manifest_sha256_and_linear_pixel_references():
    import hashlib
    manifest=json.loads((ROOT/'examples/output/visualization.json').read_text())
    assert manifest['orientation_wxyz']==[math.sqrt(.5),0.,0.,math.sqrt(.5)]
    for index,name in enumerate(['sh_original','sh_rotated']):
        record=manifest['images'][name]
        assert hashlib.sha256((ROOT/'examples/output'/record['file']).read_bytes()).hexdigest()==record['sha256']
        coefficients=manifest['irradiance_coefficients'][index]
        width,height,png=decode_png(ROOT/'examples/output'/record['file'])
        for pixel in record['reference_pixels']:
            x,y,z=pixel['normal'];a=math.sqrt(3/(4*math.pi));b=math.sqrt(15/(4*math.pi))
            basis=[1/math.sqrt(4*math.pi),-a*y,a*z,-a*x,b*x*y,-b*y*z,
                   math.sqrt(5/(16*math.pi))*(3*z*z-1),-b*x*z,math.sqrt(15/(16*math.pi))*(x*x-y*y)]
            expected=[sum(row[c]*value for row,value in zip(coefficients,basis)) for c in range(3)]
            assert pixel['irradiance_linear_rgb']==pytest.approx(expected,rel=5e-14,abs=1e-14)

            predicted=[]
            for value,albedo in zip(expected,manifest['albedo_linear_rgb']):
                light=max(0,value*albedo/math.pi)*manifest['exposure']
                mapped=light/(1+light)
                display=12.92*mapped if mapped<=.0031308 else 1.055*mapped**(1/2.4)-.055
                predicted.append(round(display*255))
            column,row=pixel['pixel'];offset=(row*width+column)*3
            assert list(png[offset:offset+3])==predicted
