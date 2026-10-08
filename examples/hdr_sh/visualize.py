"""CPU-only diffuse sphere visualization; all SH math is delegated to gem."""
import argparse
import hashlib
import json
import math
import struct
import zlib
from pathlib import Path
from gem import spherical_harmonics as sh
from gem.quaternion import Quaternion
from .environment import latlong_pixel, project_latlong
from .reference import linear_to_srgb, reflected_radiance, glsl_uniforms


def directional_environment(width=128,height=64):
    axis=[.8,0.,.6]
    image=[]
    for row in range(height):
        pixels=[]
        for column in range(width):
            direction=latlong_pixel(row,column,width,height)[2]
            peak=max(0.,sum(a*b for a,b in zip(axis,direction)))**32
            pixels.append([.08+45*peak,.06+30*peak,.04+18*peak])
        image.append(pixels)
    return image


def render_sphere(irradiance_coefficients,size=192,albedo=(.65,.65,.65),exposure=1.):
    """Orthographic camera looks along -Z; screen right +X, up +Y."""
    pixels=[];luminance=[];negative=0;references=[]
    selected={(size//2,size//2),(3*size//4,size//2),(size//2,size//4)}
    radius=.9*size/2
    for row in range(size):
        for column in range(size):
            x=(column+.5-size/2)/radius;y=(size/2-row-.5)/radius
            if x*x+y*y<=1:
                normal=[x,y,math.sqrt(max(0.,1-x*x-y*y))]
                irradiance=sh.reconstruct(irradiance_coefficients,normal)
                negative+=int(any(value<0 for value in irradiance))
                light=reflected_radiance(irradiance,albedo)
                if (column,row) in selected:
                    references.append({"pixel":[column,row],"normal":normal,
                                       "irradiance_linear_rgb":irradiance,
                                       "lambertian_linear_rgb":light})
                luminance.append((column,row,sum(a*b for a,b in zip(light,(.2126,.7152,.0722)))))
            else:light=[.015,.015,.015]
            # Display-only clipping and Reinhard tone mapping; SH data untouched.
            exposed=[max(0.,value)*exposure for value in light]
            display=[linear_to_srgb(v/(1+v)) for v in exposed]
            pixels.extend(round(min(1,max(0,v))*255) for v in display)
    weight=sum(max(0,value) for _,_,value in luminance)
    centroid=[sum(coord*max(0,value) for column,row,value in luminance
                  for coord in [column if axis==0 else row])/weight for axis in range(2)]
    stats={'size':[size,size],'minimum_byte':min(pixels),'maximum_byte':max(pixels),
           'mean_byte':sum(pixels)/len(pixels),'linear_luminance_centroid_pixels':centroid,
           'sphere_pixels':len(luminance),'negative_irradiance_pixels':negative,
           'reference_pixels':references}
    return bytes(pixels),stats


def write_png(file,width,height,pixels):
    """Minimal lossless RGB8 PNG writer; no image library or GPU required."""
    if len(pixels)!=width*height*3:raise ValueError('RGB pixel count mismatch')
    def chunk(kind,data):
        return struct.pack('>I',len(data))+kind+data+struct.pack('>I',zlib.crc32(kind+data)&0xffffffff)
    scanlines=b''.join(b'\0'+pixels[row*width*3:(row+1)*width*3] for row in range(height))
    payload=b'\x89PNG\r\n\x1a\n'+chunk(b'IHDR',struct.pack('>IIBBBBB',width,height,8,2,0,0,0))
    # Tag display output as sRGB; this is not a linear-HDR image export.
    payload+=chunk(b'sRGB',b'\0')+chunk(b'IDAT',zlib.compress(scanlines,9))+chunk(b'IEND',b'')
    Path(file).write_bytes(payload)


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path('examples/output'))
    parser.add_argument('--size',type=int,default=192)
    args=parser.parse_args(argv)
    if args.size<16:parser.error('size must be at least 16')
    args.output_dir.mkdir(parents=True,exist_ok=True)
    image=directional_environment()
    radiance=project_latlong(image)
    orientation=Quaternion([math.sqrt(.5),0.,0.,math.sqrt(.5)])
    rotated=sh.rotate_coefficients(radiance,orientation)
    coefficients=[sh.convolve_diffuse(radiance),sh.convolve_diffuse(rotated)]
    report={'schema':'gem-sh-sphere-v1','environment_size':[128,64],
            'feature_axis':[.8,0,.6],'feature_exponent':32,
            'peak_linear_rgb':[45,30,18],'ambient_linear_rgb':[.08,.06,.04],
            'randomness':'none','sphere_radius_fraction':.9,'sampling':'pixel centers',
            'rotation':'active +90 degrees about Z','orientation_wxyz':list(orientation.data),'camera':'orthographic +Z, right +X, up +Y',
            'albedo_linear_rgb':[.65]*3,'exposure':1.,'tone_mapping':'Reinhard, then sRGB',
            'radiance_coefficients':radiance,'rotated_radiance_coefficients':rotated,
            'irradiance_coefficients':coefficients,'images':{}}
    for name,coeff in zip(['sh_original','sh_rotated'],coefficients):
        pixels,stats=render_sphere(coeff,args.size)
        file=args.output_dir/(name+'.png');write_png(file,args.size,args.size,pixels)
        stats['file']=file.name;stats['file_bytes']=file.stat().st_size
        stats['sha256']=hashlib.sha256(file.read_bytes()).hexdigest()
        header=args.output_dir/(name+'_coefficients.glsl')
        header.write_text(glsl_uniforms(coeff))
        stats['coefficient_header']=header.name
        report['images'][name]=stats
    (args.output_dir/'visualization.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['images'],indent=2))


if __name__=='__main__':main()
