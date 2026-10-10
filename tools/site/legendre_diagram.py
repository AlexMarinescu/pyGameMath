"""Regenerate the documentation-only Legendre plot using the supported core API."""
from fractions import Fraction
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from gem.legendre import Legendre


def render():
    # Independent analytical checks before plotting the core evaluations.
    assert Legendre(2, 0, .5).run() == float(Fraction(-1, 8))
    assert Legendre(3, 0, .5).run() == float(Fraction(-7, 16))
    assert Legendre(4, 0, .5).run() == float(Fraction(-37, 128))
    colors = ['#253247', '#2157ab', '#ac490a', '#13754e', '#793cb5']
    rows = ['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 880 440" role="img" aria-labelledby="title desc">',
        '<title id="title">Ordinary Legendre polynomials, degrees zero through four</title>',
        '<desc id="desc">Unnormalized P_l(x) on minus one to one. All curves end at one; parity alternates by degree. Values are evaluated by gem.legendre.</desc>',
        '<rect width="880" height="440" fill="#fafbfe"/>',
        '<g font-family="sans-serif" font-size="17" fill="#253247">',
        '<text x="64" y="32" font-size="21" font-weight="600">Ordinary Legendre functions · m = 0</text>']
    for value in (-1, -.5, 0, .5, 1):
        x, y = 64+650*(value+1)/2, 62+300*(1-value)/2
        rows += [f'<path d="M{x} 62 V362 M64 {y} H714" stroke="#dbe2ec" fill="none"/>',
                 f'<text x="{x}" y="388" text-anchor="middle">{value:g}</text>',
                 f'<text x="52" y="{y+6}" text-anchor="end">{value:g}</text>']
    rows.append('<path d="M389 62 V362 M64 212 H714" stroke="#8493a9" fill="none"/>')
    for degree, color in enumerate(colors):
        coords = []
        for i in range(401):
            argument = -1+2*i/400
            value = Legendre(degree, 0, argument).run()
            coords.append(f'{64+650*(argument+1)/2:.3f},{62+300*(1-value)/2:.3f}')
        dash = ['', '10 3', '3 3', '10 3 3 3', '14 5'][degree]
        rows += [f'<polyline points="{" ".join(coords)}" fill="none" stroke="{color}" stroke-width="2.7" stroke-dasharray="{dash}"/>',
                 f'<path d="M746 {85+degree*38} H785" stroke="{color}" stroke-width="3" stroke-dasharray="{dash}"/>',
                 f'<text x="796" y="{91+degree*38}">P{degree}</text>']
    rows += ['<text x="389" y="422" text-anchor="middle">x · independent variable</text>', '</g></svg>']
    return '\n'.join(rows)+'\n'


if __name__ == '__main__':
    target = ROOT/'docs/assets/diagrams/legendre.svg'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(render())
    print(target)
