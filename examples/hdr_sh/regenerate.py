"""Regenerate all committed example images and machine-readable references."""
import argparse
from pathlib import Path
from . import reference, visualize


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,default=Path('examples/output'))
    args=parser.parse_args();args.output_dir.mkdir(parents=True,exist_ok=True)
    reference.main(['--output',str(args.output_dir/'reference.json'),
                    '--glsl-output',str(args.output_dir/'irradiance_coefficients.glsl')])
    visualize.main(['--output-dir',str(args.output_dir)])



if __name__=='__main__':main()
