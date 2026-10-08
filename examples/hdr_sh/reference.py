"""Run with python -m examples.hdr_sh.reference; no graphics context required."""
import argparse
import json
import math
from pathlib import Path
from gem import spherical_harmonics as sh
from gem.quaternion import Quaternion
from .environment import load_fixture, load_raw, load_rgbe, project_latlong


def linear_to_srgb(value):
    """Final display encoding only; radiance/irradiance computations stay linear."""
    return 12.92*value if value <= .0031308 else 1.055*value**(1/2.4)-.055


def reflected_radiance(irradiance, albedo):
    return [light*reflectance/math.pi for light,reflectance in zip(irradiance,albedo)]


def run(image, mapping, orientation, albedo):
    if len(albedo)!=3 or not all(math.isfinite(v) and 0<=v<=1 for v in albedo):
        raise ValueError("albedo must be three finite values within [0,1]")
    if mapping=='angular':coefficients=sh.project_angular_probe(image)
    elif mapping=='latlong':coefficients=project_latlong(image)
    else:raise ValueError("mapping must be angular or latlong")
    rotated=sh.rotate_coefficients(coefficients,orientation)
    irradiance=sh.convolve_diffuse(rotated)
    directions=[[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[-1.,0.,0.],[0.,-1.,0.],[0.,0.,-1.],[.3,.4,math.sqrt(.75)]]
    evaluations=[]
    for direction in directions:
        light=sh.reconstruct(irradiance,direction)
        reflected=reflected_radiance(light,albedo)
        evaluations.append({'normal':direction,'irradiance_linear_rgb':light,
                            'lambertian_linear_rgb':reflected,
                            'display_srgb': [linear_to_srgb(v) for v in reflected]})
    return {'schema':'gem-sh-reference-v1','basis':'canonical real SH, Condon-Shortley',
            'index':'l*(l+1)+m','layout':'coefficient-major RGB; nine vec3 rows',
            'mapping':mapping,'image_size':[len(image[0]),len(image)],
            'working_color':'linear RGB; display encoding assumes linear sRGB primaries',
            'orientation_wxyz':list(orientation.data),'rotation':'active: f_rotated(d)=f_original(R^-1 d)',
            'albedo_linear_rgb':list(albedo),'radiance_coefficients':coefficients,
            'rotated_radiance_coefficients':rotated,'irradiance_coefficients':irradiance,
            'evaluations':evaluations}


def glsl_uniforms(coefficients):
    rows=["    vec3("+", ".join(format(v,'.9e') for v in row)+")" for row in coefficients]
    return 'const vec3 referenceIrradianceSH[9] = vec3[9](\n'+',\n'.join(rows)+'\n);\n'


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path(__file__).parent/'fixtures/asymmetric.json')
    parser.add_argument('--format',choices=['fixture','raw','rgbe'],default='fixture')
    parser.add_argument('--mapping',choices=['angular','latlong'])
    parser.add_argument('--width',type=int);parser.add_argument('--height',type=int)
    parser.add_argument('--quaternion',type=float,nargs=4,default=[math.sqrt(.5),0,0,math.sqrt(.5)],metavar=('W','X','Y','Z'))
    parser.add_argument('--albedo',type=float,nargs=3,default=[.5,.5,.5])
    parser.add_argument('--output',type=Path)
    parser.add_argument('--glsl-output',type=Path)
    args=parser.parse_args(argv)
    if args.format=='fixture':
        if args.mapping not in (None,'latlong'):parser.error('synthetic fixture is latitude-longitude')
        args.mapping='latlong'
        image=load_fixture(args.input)
    elif args.format=='raw':
        if args.mapping is None:parser.error('real input requires explicit --mapping')
        if args.width is None or args.height is None:parser.error('raw input needs --width and --height')
        image=load_raw(args.input,args.width,args.height)
    else:
        if args.mapping is None:parser.error('real input requires explicit --mapping')
        image=load_rgbe(args.input)
    result=run(image,args.mapping,Quaternion(args.quaternion),args.albedo)
    text=json.dumps(result,indent=2,allow_nan=False)+'\n'
    if args.output:args.output.write_text(text)
    else:print(text,end='')
    if args.glsl_output:args.glsl_output.write_text(glsl_uniforms(result['irradiance_coefficients']))


if __name__=='__main__':main()
