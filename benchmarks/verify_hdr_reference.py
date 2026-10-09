"""Regenerate HDR/SH outputs without overwriting committed references."""
import argparse
import hashlib
import json
from pathlib import Path
import struct
import subprocess
import sys
import zlib

ROOT=Path(__file__).resolve().parents[1]


def png_rgb(path):
    data=path.read_bytes();assert data[:8]==b'\x89PNG\r\n\x1a\n'
    position=8;compressed=b''
    while position<len(data):
        length,=struct.unpack('>I',data[position:position+4]);kind=data[position+4:position+8]
        content=data[position+8:position+8+length];position+=12+length
        if kind==b'IHDR':width,height,depth,color,*_=struct.unpack('>IIBBBBB',content);assert (depth,color)==(8,2)
        if kind==b'IDAT':compressed+=content
    raw=zlib.decompress(compressed);stride=width*3+1
    assert len(raw)==height*stride and all(raw[row*stride]==0 for row in range(height))
    return width,height,b''.join(raw[row*stride+1:(row+1)*stride] for row in range(height))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args();args.output_dir=args.output_dir.resolve()
    assert args.output_dir!=(ROOT/'examples/output').resolve(), 'Use a separate regeneration directory'
    command=[sys.executable,'-m','examples.hdr_sh.regenerate','--output-dir',str(args.output_dir)]
    subprocess.run(command,cwd=ROOT,check=True,capture_output=True,text=True)
    records={}
    for original in sorted((ROOT/'examples/output').iterdir()):
        if not original.is_file():continue
        generated=args.output_dir/original.name
        assert generated.exists(),original.name
        old,new=original.read_bytes(),generated.read_bytes()
        row={'committed_sha256':hashlib.sha256(old).hexdigest(),'generated_sha256':hashlib.sha256(new).hexdigest(),
             'bytes_identical':old==new,'bytes':len(new)}
        assert old==new, original.name
        if original.suffix=='.png':
            w,h,pixels=png_rgb(generated);ow,oh,reference=png_rgb(original)
            assert (w,h,pixels)==(ow,oh,reference)
            row.update(width=w,height=h,decoded_rgb_identical=True,minimum_byte=min(pixels),maximum_byte=max(pixels),
                       mean_byte=sum(pixels)/len(pixels),decoded_rgb_sha256=hashlib.sha256(pixels).hexdigest())
        records[original.name]=row
    old=json.loads((ROOT/'examples/output/visualization.json').read_text())
    new=json.loads((args.output_dir/'visualization.json').read_text())
    fields=('radiance_coefficients','rotated_radiance_coefficients','irradiance_coefficients')
    equal={field:old[field]==new[field] for field in fields};assert all(equal.values())
    for name in ('sh_original','sh_rotated'):
        assert old['images'][name]['reference_pixels']==new['images'][name]['reference_pixels']
    result={'command':command,'files':records,'coefficient_equality':equal,'selected_linear_pixels_identical':True,
            'goldens_modified':False,'environment':{'python':sys.version,'zlib':zlib.ZLIB_VERSION,'zlib_runtime':zlib.ZLIB_RUNTIME_VERSION},
            'limitations':'byte identity verified on this host; Python/libm and zlib/platform changes can affect floating outputs or compressed bytes'}
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(f'{len(records)} regenerated files byte-identical; coefficients, selected linear values and decoded pixels match')


if __name__=='__main__':main()
