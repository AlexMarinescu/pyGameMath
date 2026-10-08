"""Example-only image adapters; decoded values remain linear RGB."""
import json
import math
import struct
from pathlib import Path
from gem import spherical_harmonics as sh
from gem.vector import Vector


def latlong_pixel(row, column, width, height):
    theta = math.pi*(row+0.5)/height
    phi = 2*math.pi*(column+0.5)/width
    direction = [math.sin(theta)*math.cos(phi), math.sin(theta)*math.sin(phi), math.cos(theta)]
    weight = 2*math.pi/width*(math.cos(math.pi*row/height)-math.cos(math.pi*(row+1)/height))
    return theta, phi, direction, weight


def project_latlong(image):
    """Project center radiance with exact pixel-area weights, not angular weights."""
    if not image or not image[0] or any(len(row) != len(image[0]) for row in image):
        raise ValueError("require a nonempty rectangular RGB image")
    height, width = len(image), len(image[0])
    samples, radiances, weights = [], [], []
    for row in range(height):
        for column in range(width):
            theta, phi, direction, weight = latlong_pixel(row,column,width,height)
            sample = sh.SPHSample(theta,phi,Vector(3,direction),9)
            sample.values = [sh.SPH(l,m,theta,phi) for l in range(3) for m in range(-l,l+1)]
            samples.append(sample);radiances.append(image[row][column]);weights.append(weight)
    return sh.project_radiance(samples,radiances,weights)


def synthetic_environment(width, height, constant, linear_xyz):
    """Linear RGB function L(d)=constant+linear_xyz*d sampled at pixel centers."""
    if width <= 0 or height <= 0:
        raise ValueError("positive image dimensions required")
    return [[[constant[c]+sum(a*b for a,b in zip(linear_xyz[c],latlong_pixel(row,col,width,height)[2]))
              for c in range(3)] for col in range(width)] for row in range(height)]


def load_fixture(file):
    data=json.loads(Path(file).read_text())
    if data['format'] != 'synthetic-linear-latlong':
        raise ValueError("unsupported synthetic fixture format")
    return synthetic_environment(data['width'],data['height'],data['constant'],data['linear_xyz'])


def load_raw(file, width, height):
    if width <= 0 or height <= 0:
        raise ValueError("positive raw image dimensions required")
    count=width*height*3
    with open(file,'rb') as handle:data=handle.read(count*4)
    if len(data) != count*4:raise ValueError("truncated raw RGB float32 image")
    values=struct.unpack('%df'%count,data)
    return [[list(values[(r*width+c)*3:(r*width+c+1)*3]) for c in range(width)] for r in range(height)]


def load_rgbe(file):
    """Read flat/new scanline-RLE RGBE with standard -Y/+X orientation.

    This intentionally small reader excludes XYZE, legacy run markers and
    other scan orders. Stored RGB is returned without display/exposure changes.
    """
    with open(file,'rb') as handle:
        def read(count):
            result=handle.read(count)
            if len(result) != count:raise ValueError("truncated RGBE file")
            return result
        if handle.readline().strip() not in (b'#?RADIANCE',b'#?RGBE'):
            raise ValueError("not a Radiance RGBE file")
        header=[]
        while True:
            line=handle.readline()
            if not line:raise ValueError("truncated RGBE header")
            if not line.strip():break
            header.append(line.strip())
        if b'FORMAT=32-bit_rle_rgbe' not in header:
            raise ValueError("only RGBE encoding is supported")
        resolution=handle.readline().split()
        if len(resolution)!=4 or resolution[0]!=b'-Y' or resolution[2]!=b'+X':
            raise ValueError("RGBE reader requires -Y height +X width scan order")
        height,width=int(resolution[1]),int(resolution[3])
        if min(width,height)<=0:raise ValueError("positive RGBE dimensions required")
        image=[]
        for _ in range(height):
            first=read(4)
            if 8 <= width <= 32767 and first[:2]==b'\x02\x02' and first[2]<128:
                if (first[2]<<8)+first[3]!=width:raise ValueError("RGBE scanline width mismatch")
                channels=[]
                for _ in range(4):
                    channel=[]
                    while len(channel)<width:
                        code=read(1)[0]
                        if code==0:raise ValueError("zero-length RGBE run")
                        count=code-128 if code>128 else code
                        if len(channel)+count>width:raise ValueError("RGBE run exceeds scanline")
                        channel.extend([read(1)[0]]*count if code>128 else read(count))
                    channels.append(channel)
                pixels=list(zip(*channels))
            else:
                data=first+read(4*(width-1))
                pixels=[data[i:i+4] for i in range(0,len(data),4)]
                if any(p[:3]==b'\x01\x01\x01' for p in pixels):
                    raise ValueError("legacy RGBE run markers are unsupported")
            row=[]
            for red,green,blue,exponent in pixels:
                scale=math.ldexp(1.,exponent-136) if exponent else 0.
                row.append([red*scale,green*scale,blue*scale])
            image.append(row)
    return image
