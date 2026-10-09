"""Recheck selected cases over interleaved rounds after a noisy full-suite run."""
import argparse
import json
from pathlib import Path

from core_baseline import prepare, measure


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    calls, _ = prepare()
    names = ['vector2_dot', 'matrix3_raw_inverse', 'matrix4_raw_inverse',
             'sh_project_n1024_b3', 'sh_rotate_l2_scalar', 'ray_duplicate',
             'ray_translate', 'ray_matrix_rotate', 'ray_quaternion_rotate']
    rounds = []
    for _ in range(3):
        rounds.append({name: measure(calls[name], 7, .05) for name in names})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({'methodology': 'three interleaved rounds; 7 trials, .05s calibration per case',
                                     'rounds': rounds}, indent=2) + '\n')
    print(f'{len(names)} cases, 3 rounds saved to {args.output}')


if __name__ == '__main__':
    main()
