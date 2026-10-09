# pyGameMath · gem

**Pure-Python graphics and computational mathematics.**

gem is a lightweight mathematics library for graphics, geometry, game development,
simulation and computational research. Use it for small vectors and matrices,
orientation, curves and CPU-side lighting references, independently of an engine
or renderer. It is not a game engine or an array-oriented replacement for NumPy.

Mathematical algorithms are implemented in Python, with standard-library `math`
and `ctypes` interoperability. There is no NumPy, Cython or compiled numerical
backend requirement. The existing runtime dependency is **six**.

## Project status

The repository is undergoing [1.0 modernization](ROADMAP.md). Current master
contains audited correctness fixes, numerical-stability improvements, focused
performance work and supported core modules. Documentation and release engineering
are still in progress; a 1.0 release and blanket API stability are not claimed.

The examples below target **current repository code**, not the historical PyPI
release. CPython 3.12 on Linux is the verified reference environment. Other
interpreters/platforms, including the deliberate Python 2.7 legacy target, need
explicit assessment and testing before support is advertised.

## Implemented features

| Area | Available today |
|---|---|
| Algebra | Vectors, dot/cross products, stable lengths/normalization, Hamilton quaternions |
| Matrices and transforms | Small-matrix arithmetic, 2x2/3x3/4x4 determinants/inverses, rotation, translation, shear, camera and project/unproject helpers |
| Geometry | Plane construction/normalization, barycentric weights, ray duplication and rigid transforms; no ray intersection API |
| Curves | Scalar/Vector quadratic and cubic Bezier evaluation, cubic paths and depth-limited adaptive sampling |
| Polynomial basis | Ordinary/associated Legendre functions with documented phase and domains |
| Lighting | Real SH basis, RGB radiance projection, diffuse convolution/reconstruction and analytical scalar/RGB rotation through L2 |

Use canonical `gem.bezier`, `gem.legendre` and `gem.spherical_harmonics` imports.
Historical experimental paths remain thin [compatibility reexports](docs/EXPERIMENTAL_MIGRATION.md).

## Installation

The distribution name and import namespace are both **`gem`**. To use the
modernized development code, install this repository in a fresh environment:

```sh
git clone https://github.com/AlexMarinescu/pyGameMath.git
cd pyGameMath
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install .
python -c "from gem.vector import Vector; print((Vector(3, [1, 2, 3]) + Vector(3, [4, 5, 6])).vector)"
```

The import check prints `[5, 7, 9]`. These shell commands use POSIX activation;
see [installation](docs/getting-started/installation.md) for Windows commands,
development edits and isolated verification details.

[PyPI](https://pypi.org/project/gem/) lists the historical 0.1.12 release from
2017. Installing it is not the installation path for these current features.
Repository metadata still declares the same version; use a fresh environment
and the source URL/check-out to identify the code. Packaging/version modernization
is separate work, and no new release has been published here.

## Quick start

### Vectors

```python
from gem.vector import Vector, cross

a = Vector(3, [3.0, 4.0, 0.0])
print(a.magnitude())                  # 5.0
print(a.dot(Vector(3, [1, 0, 0])))    # 3.0
print(cross(Vector(3, [1, 0, 0]), Vector(3, [0, 1, 0])).vector)
# [0, 0, 1]
assert a.normalize().vector == [0.6, 0.8, 0.0]
assert a.vector == [3.0, 4.0, 0.0]    # returning normalization preserves input
```

### Quaternion rotation

```python
from gem.quaternion import quat_from_axis_angle, quat_rotate_vector
from gem.vector import Vector

rotation = quat_from_axis_angle([0, 0, 1], 90)  # degrees; [w,x,y,z] storage
turned = quat_rotate_vector(rotation, Vector(3, [1, 0, 0]))
print([round(value, 6) for value in turned.vector])  # [0.0, 1.0, 0.0]
assert abs(turned.vector[1] - 1.0) < 1e-14
```

These are dimensioned `Vector`/`Matrix` classes, not separate Vector3/Matrix4
constructors. Matrices use row-vector mathematics despite the `Matrix * Vector`
syntax, with translation in the final row. Angle units vary by API. The
[practical quick start](docs/getting-started/quick-start.md) covers matrix
composition/inversion, Bezier curves, SH and ctypes with executable examples;
[conventions](docs/architecture/conventions.md) explains the historical distinctions.

## Practical uses

Use camera/object transforms and quaternion orientation for graphics integration;
Bezier paths for procedural curves; and the SH core as an independent CPU
reference for diffuse environment-lighting calculations. Matrix float32 ctypes
snapshots can interface with foreign APIs, provided their layout/transpose
requirements are handled explicitly.

The [headless HDR/SH example](examples/hdr_sh/README.md) generates an asymmetric
linear-HDR environment in code, projects canonical L2 RGB coefficients, rotates
lighting analytically, convolves once and renders CPU diffuse-sphere PNGs. It also
exports coefficients and GLSL reference expressions. It needs no GPU/OpenGL
context and does not claim GPU execution.

| Original lighting | Active +90-degree Z lighting rotation |
|---|---|
| ![CPU diffuse sphere with original lighting](examples/output/sh_original.png) | ![CPU diffuse sphere with rotated lighting](examples/output/sh_rotated.png) |

Regenerate into a separate directory from the checkout:

```sh
python -m examples.hdr_sh.regenerate --output-dir /tmp/gem-hdr-reference
```

The example source is not installed as a core renderer. Coefficients and selected
linear pixel values provide numerical evidence in addition to the images.

## Correctness and performance

Independent known answers, algebraic/geometric invariants, extreme-scale cases
and ownership tests define supported behavior. Stable lengths and scaled Matrix3/4
inverses avoid avoidable overflow/underflow, but do not guarantee every operation
for arbitrary nonfinite or ill-conditioned inputs. Approximation and integration
limits are explicit; Vector equality remains exact.

[Performance monitoring](benchmarks/REGRESSION_POLICY.md) retains deterministic
workloads, repeated timings and environment metadata. Comparisons flag possible
regressions for review rather than impose brittle timing gates on shared CI.
Pure Python and numerical contracts take precedence over unsupported speed claims.

## Documentation

| Start here | Deeper reference |
|---|---|
| [Getting started](docs/getting-started/README.md) | [API reference](docs/api/index.md) |
| [Installation](docs/getting-started/installation.md) | [Mathematical conventions](docs/architecture/conventions.md) |
| [Practical tutorials](docs/tutorials/index.md) | [Architecture charter](docs/architecture/philosophy.md) |
| [Documentation index](docs/README.md) | [Compatibility policy](docs/architecture/compatibility.md) |
| [HDR/SH reference](examples/hdr_sh/README.md) | [Benchmark workflow](benchmarks/README.md) |

The topic reference documents current signatures, domains and ownership, with
[declaration coverage](docs/api/coverage.md). No generated documentation site exists.
Historical Wiki content has not yet been migrated.

## Roadmap

[ROADMAP.md](ROADMAP.md) is the single source of truth. Proposed future work
includes geometry foundations, procedural noise/sampling, SDFs and spatial
mathematics, volumetric/lighting mathematics and advanced computational geometry.
Those expansions are **planned**, not shipped; later version assignments are
provisional with no release dates.

## Compatibility and dependencies

Retain pure-Python operation, standard-library integrations and six. Current
verification covers CPython 3.12.14/Linux x86_64; older classifiers are not verified
support claims. Python 2.7 remains a legacy engineering target with known source/
tooling blockers, not tested compatibility. See the
[policy and evidence](docs/architecture/compatibility.md) for numeric reproducibility,
interpreter assessment and migration obligations.

Constructors may retain caller-owned lists, while returning operations normally
allocate fresh results. Read ownership and unit prerequisites before integrating
with existing mutable application state.

## Contributing

Start with the [roadmap and open decisions](docs/development/roadmap.md), then
add independent mathematical references and preserve established APIs, units and
ownership. In an activated development environment:

```sh
python -m pip install -r requirements-audit.txt
python -m pytest -q
python tools/check_architecture_docs.py --examples
```

The pinned audit tools require modern Python. See
[documentation verification](docs/getting-started/verification.md) and the
[benchmark policy](benchmarks/REGRESSION_POLICY.md) for reproducible checks.

## License and attribution

Created by Alex Marinescu, originally for learning graphics mathematics and
personal OpenGL projects. Distributed under the [BSD 2-Clause license](LICENSE),
copyright 2015–2026 Alex Marinescu. Historical attribution remains intact.

Explore the [headless graphics gallery](docs/examples/index.md) for reproducible visual examples.
