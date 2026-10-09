# Executed gallery verification

The branch starts at merged PR #42, master
`da2ad2f7715c81851d7853031616036ddfb4f9e5`. Verification used CPython 3.12.14,
Linux x86_64, six 1.17.0 and pytest 9.1.1. Build tooling was setuptools 84.0.0,
wheel 0.48.0 and packaging 26.3; none is a new gem runtime dependency.

## Results

| Check | Result |
| --- | --- |
| Full `python -m pytest -q` | 2,257 passed; 0 failed, 0 xfailed, 0 skipped (9.96s) |
| New gallery regressions | 8 passed, including corruption detection and repeatability |
| Phase 4A–4D checks plus gallery links | 42 Python blocks executed; 268 declarations and 7 constants verified |
| Source, wheel and sdist documentation checks | Passed; installed module paths asserted outside the checkout |
| Protected scope | All 93 pre-existing mathematical/test/example/benchmark and packaging/license files unchanged |
| Clean checkout regeneration | Four SVGs and JSON match committed bytes |
| Isolated wheel and sdist regeneration | All five artifacts match committed bytes |
| Existing HDR regeneration | Both PNG SHA-256 hashes match; golden files unchanged |
| Visual inspection | All four SVGs rasterized for review; labels and shapes inspected; existing HDR pair inspected |

The [machine-readable results](verification-results.json) record measurements
and checks. The [numerical manifest](../../examples/showcase/output/measurements.json)
contains the SVG hashes and geometric values.

| Output | Format / dimensions | Bytes |
| --- | --- | ---: |
| `examples/showcase/output/vectors.svg` | SVG, 960×560 | 3,302 |
| `examples/showcase/output/transforms.svg` | SVG, 960×560 | 5,803 |
| `examples/showcase/output/quaternions.svg` | SVG, 960×560 | 4,774 |
| `examples/showcase/output/bezier.svg` | SVG, 960×560 | 9,235 |
| `examples/showcase/output/measurements.json` | UTF-8 JSON | 8,967 |
| `examples/output/sh_original.png` (existing) | RGB8 sRGB PNG, 192×192 | 11,407 |
| `examples/output/sh_rotated.png` (existing) | RGB8 sRGB PNG, 192×192 | 16,733 |

Independent checks establish sum=(2,3,0), dot=-1, cross=(0,0,7); six transformation
coordinate mappings; five SLERP angles; quadratic/cubic midpoints; curve endpoints
and monotonic parameter order. Cubic adaptive paths contain 5 and 14 points.
Against 1,001 independent expanded-polynomial reference values, maximum distances
are 0.2504016663 and 0.0299070097, below the illustrative 0.4 and 0.05 distances.
These are measured example results, not universal subdivision guarantees.

## Commands

From the checkout:

```sh
python -m pytest -q
python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --output /tmp/phase4e-docs.json
python -m examples.showcase.regenerate
python -m examples.showcase.regenerate --verify
```

Clean source verification copied only the branch changes into a local clone of
the merged master, then ran these commands under `/tmp/phase4e-source`:

```sh
/tmp/phase4e-wheel-env/bin/python setup.py sdist bdist_wheel
/tmp/phase4e-wheel-env/bin/python -m pip install --no-index dist/gem-0.1.12-py3-none-any.whl
/tmp/phase4e-wheel-env/bin/python -m examples.showcase.regenerate --verify
/tmp/phase4e-wheel-env/bin/python -m examples.showcase.regenerate --output-dir /tmp/phase4e-clean-gallery
/tmp/phase4e-wheel-env/bin/python -m examples.hdr_sh.regenerate --output-dir /tmp/phase4e-clean-hdr
```

The wheel and sdist were installed into separate fresh venvs with six. Example
folders were copied into `/tmp/phase4e-wheel-examples`, without gem. From `/tmp`,
isolated wheel execution used:

```sh
/tmp/phase4e-wheel-env/bin/python -I -c 'import gem,sys; assert gem.__file__.startswith("/tmp/phase4e-wheel-env/lib/python3.12/site-packages/"); sys.path.insert(0,"/tmp/phase4e-wheel-examples"); from showcase.regenerate import main; main(["--output-dir","/tmp/phase4e-installed-gallery"])'
/tmp/phase4e-wheel-env/bin/python -I -c 'import gem,sys; assert gem.__file__.startswith("/tmp/phase4e-wheel-env/lib/python3.12/site-packages/"); sys.path.insert(0,"/tmp/phase4e-wheel-examples"); from hdr_sh.regenerate import main; sys.argv=["regenerate","--output-dir","/tmp/phase4e-installed-hdr"]; main()'
/tmp/phase4e-wheel-env/bin/python -I /workspace/pyGameMath/tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --package-root /tmp/phase4e-wheel-env/lib/python3.12/site-packages --output /tmp/phase4e-wheel-docs.json
```

The equivalent sdist venv ran gallery regeneration to `/tmp/phase4e-sdist-gallery`
and the same documentation checker with its own installed package root. All five
files were compared byte-for-byte in all three contexts. Existing PNG hashes:

- Original: `dfc2e0f5aa7cdee17dc32231ae3258c551cc886eff2bbb01cc43f2e38360174a`
- Rotated: `c9b205dc7dee31d731685d996326badf9da62c48fe69713b6fbe2fea4fccaac7`

SVG review used local Inkscape rasterization into `/tmp/showcase-review`; it is
inspection tooling only. Regeneration requires no Inkscape, Pillow, image decoder,
GPU, downloaded asset or external font. SVG font appearance and floating-point
JSON serialization can vary across platforms; independent numerical and structural
checks complement hashes. No display-image bytes are the sole correctness oracle.

## Compatibility and open decisions

There are no core, API, mathematical, dependency, release-metadata or license
changes. Existing experimental compatibility paths remain intact. The examples
are source-checkout tools rather than new installed runtime modules. No roadmap
functionality or documentation website is implemented. No new mathematical
contract decision is needed. Future website styling and font selection remain
presentation work; the [roadmap](../../ROADMAP.md) remains canonical.

[Gallery](index.md) · [Asset manifest](ASSET-MANIFEST.md)
