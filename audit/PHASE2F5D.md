# HDR-to-SH shader reference

The complete example lives under `examples/hdr_sh`, leaving gem and its
mandatory dependencies unchanged. It supports native-endian float32 angular
probes, a separately defined latitude-longitude adapter and a limited pure-
Python Radiance RGBE decoder. The canonical core APIs handle basis evaluation,
projection, analytical rotation, cosine convolution and reconstruction.

The primary numeric pipeline exports nine canonical RGB rows, known-normal
comparisons and GLSL constants. The CPU visualizer generates a sharp asymmetric
HDR environment in code, projects once and renders the same diffuse sphere
before/after a +90° Z coefficient rotation. Pixel values stay linear until
exposure/tone mapping/display encoding. The two actual PNGs and their manifest
are committed. No GPU context, image service, downloaded image or manually
constructed output is involved.

See [formats, mappings, shader use and limitations](../examples/hdr_sh/README.md).
The Radiance reader's format/scan-order limits are explicit; angular disks,
latitude-longitude maps and mirrored-ball photographs are not conflated.
Legacy coefficient arrays are never silently treated as canonical.

## Executed commands

From `/workspace/pyGameMath`:

```sh
/workspace/.venvs/pyGameMath/bin/python -m examples.hdr_sh.regenerate
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/hdr-example-final.xml
```

A core wheel was built from a clean temporary copy of gem/setup.py and pip
installed with --no-index --no-deps. Smoke execution explicitly checked that
gem.spherical_harmonics resolved under `/tmp/hdr-wheel-installed`, not source.
A fresh virtual environment `/tmp/hdr-clean-venv` then installed that wheel.
The existing mandatory six 1.17.0 pure-Python module and distribution metadata
were copied offline into that environment; no image dependencies were present.
Only example source/fixtures were copied to `/tmp/hdr-clean/examples`, with no
reference outputs or core source there before the first regeneration.

From `/tmp/hdr-clean` the following commands were executed:

```sh
/tmp/hdr-clean-venv/bin/python -m examples.hdr_sh.regenerate
sha256sum examples/output/sh_original.png examples/output/sh_rotated.png
cmp examples/output/sh_original.png /workspace/pyGameMath/examples/output/sh_original.png
cmp examples/output/sh_rotated.png /workspace/pyGameMath/examples/output/sh_rotated.png
cmp examples/output/visualization.json /workspace/pyGameMath/examples/output/visualization.json
cmp examples/output/reference.json /workspace/pyGameMath/examples/output/reference.json
cmp examples/output/sh_original_coefficients.glsl /workspace/pyGameMath/examples/output/sh_original_coefficients.glsl
cmp examples/output/sh_rotated_coefficients.glsl /workspace/pyGameMath/examples/output/sh_rotated_coefficients.glsl
```

All comparisons pass. Both hashes match the committed images:

| Image | SHA-256 | Bytes |
|---|---|---:|
| sh_original.png | dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a | 11407 |
| sh_rotated.png | c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7 | 16733 |

Both images are 192x192 RGB8 with byte range 0–211 and mean 102.34946469907408.
The positive linear-luminance centroid moves from (122.8861,95.5) to
(95.5,68.1139), correctly moving illumination from right to up. Each image has
23436 sphere pixels and 300 negative-irradiance pixels from L2 ringing before
display clipping. PNGs were also visually inspected and CRC/pixel data checked.
Python 3.12.14, six 1.17.0 and zlib 1.3.2 are used. Cross-platform libm and
compression differences can affect byte-identical hashes; raw numerical
comparisons retain ordinary precision tolerances.

## Regression results

Base: 6233355f50f5dc05618026d4b533fe53b977cc13, merged PR #29.
Baseline: 1752 passed, 5 xfailed. Final: **1771 passed, 0 failed, 5 xfailed**.
Exact unrelated expected-failure identities are unchanged. The 19 new cases
cover exact latlong pixel areas/orientation, constant/asymmetric/directional
lighting, known rotations, channel independence, shader float32 expressions,
Lambertian albedo/pi, linear/display behavior, RGBE decoding/errors, raw angular
integration, ownership, controlled L0/L1/rotated L1 rendering and committed
pixel/manifest checks. Independent Cartesian polynomials validate linear pixel
values before tone mapping; PNGs are not the sole oracle.

No core implementation, renderer, GI, visibility or experimental algorithm is
changed. Examples are not installed as gem runtime packages; copy/run their
source separately when validating an installed wheel.
