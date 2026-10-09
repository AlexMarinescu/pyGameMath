"""Rebuild deterministic SVG showcase diagrams and numerical evidence offline."""
import argparse
import hashlib
import json
from pathlib import Path
from .scenes import SCENES
from .verify import check_numerics, check_svg


def generate():
    artifacts, data, metadata = {}, {}, {}
    for name, scene in SCENES.items():
        svg, data[name] = scene()
        dimensions = check_svg(svg)
        payload = svg.encode('utf-8')
        artifacts[name+'.svg'] = payload
        metadata[name] = dict(file=name+'.svg', dimensions=dimensions, bytes=len(payload),
                              sha256=hashlib.sha256(payload).hexdigest())
    numerical = check_numerics(data)
    report = dict(schema='gem-showcase-v1', randomness='none', color_space='SVG CSS sRGB',
                  artifacts=metadata, scenes=data, independent_verification=numerical)
    artifacts['measurements.json'] = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()
    return artifacts


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path(__file__).with_name('output'))
    parser.add_argument('--verify', action='store_true', help='compare without writing files')
    args = parser.parse_args(argv)
    artifacts = generate()
    if not args.verify:
        args.output_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in artifacts.items():
        destination = args.output_dir/name
        if args.verify:
            if destination.read_bytes() != payload:
                raise SystemExit('artifact mismatch: '+str(destination))
        else:
            destination.write_bytes(payload)
        print(('verified ' if args.verify else 'wrote ')+str(destination)+' ('+str(len(payload))+' bytes)')


if __name__ == '__main__':
    main()
