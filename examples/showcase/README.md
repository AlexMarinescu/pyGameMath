# Headless graphics showcase

From the repository root, with gem's existing `six` dependency installed:

```sh
python -m examples.showcase.regenerate
python -m examples.showcase.regenerate --verify
```

No GPU, window server, downloads or image packages are needed. The command writes
only `vectors.svg`, `transforms.svg`, `quaternions.svg`, `bezier.svg` and
`measurements.json` in `examples/showcase/output/`. Use `--output-dir /tmp/gem-gallery`
to generate elsewhere. Verification rebuilds in memory and compares without writes.

`scenes.py` uses gem for vector geometry, matrices, quaternion SLERP and rotation,
and Bezier evaluation/subdivision. `svg.py` only maps drawing coordinates and
writes escaped SVG elements. `verify.py` checks independent coordinate formulas,
known angles, endpoints and dense curve-reference distances before files are written.
It is example validation, not another implementation of the core algorithms.

The SVG viewBox is 960×560, with CSS sRGB colors and generic sans-serif labels.
Geometry serializes to six decimal places. Viewer font selection can change text
appearance. SHA-256 hashes describe text artifacts, not identical rasterization
across platforms. Python/libm rounding can affect the full-precision JSON values;
numerical checks remain the primary correctness test.

The existing HDR example is reused directly; this command never touches
`examples/output/`. To verify that reference separately:

```sh
python -m examples.hdr_sh.regenerate --output-dir /tmp/gem-gallery-hdr
```

See the [gallery](../../docs/examples/index.md) and
[asset manifest](../../docs/examples/ASSET-MANIFEST.md). Examples stay outside the
installed runtime package. To exercise a wheel, install it and six into a clean
venv, copy `examples/showcase/` into `/tmp/gem-gallery-wheel/showcase/`, then run
`python -I -c 'import sys; sys.path.insert(0,"/tmp/gem-gallery-wheel"); from showcase.regenerate import main; main(["--output-dir","/tmp/gem-gallery-wheel/output"])'`
from `/tmp` with that venv's Python. Only the examples directory is added to the
search path; gem resolves to the installed wheel, not the source tree.
