<div class="gem-hero" markdown>

<p class="gem-kicker">pyGameMath · the gem package</p>

# Graphics mathematics, in Python

Vectors, transformations, rotations, curves and lighting—small, readable tools
for understanding and using graphics mathematics. The `gem` package is pure
Python, with standard-library ctypes integration and six as its only runtime
dependency.

<div class="gem-actions" markdown>

[Get started](getting-started/quick-start.md){ .md-button .md-button--primary }
[Explore the tutorials](tutorials/index.md){ .md-button }
[Browse the API](api/index.md){ .md-button }

</div>

</div>

## Your first calculation

Install the current repository to use the audited development code:

```sh
python -m pip install "git+https://github.com/AlexMarinescu/pyGameMath.git"
```

```python
from gem.vector import Vector

velocity = Vector(3, [3, 4, 0])
direction = velocity.normalize()  # A fresh unit direction; velocity is preserved.
assert velocity.magnitude() == 5.0
assert direction.vector == [0.6, 0.8, 0.0]
assert velocity.vector == [3, 4, 0]
```

Follow the [installation guide](getting-started/installation.md) for a virtual
environment and compatibility details, then try the [quick start](getting-started/quick-start.md).

## A practical mathematical core

<div class="gem-cards" markdown>

<div class="gem-card" markdown>

### Vectors and algebra

Component arithmetic, dot and cross products, stable norms, matrix products
and inverses. Start with displacement and direction.

[Vectors tutorial](tutorials/vectors.md) · [Vector API](api/vector.md) · [Matrix API](api/matrix.md)

</div>
<div class="gem-card" markdown>

### Transformations and rotations

Compose row-vector transforms, follow camera coordinates and interpolate
quaternion orientations. See the order of operations in the results.

[Transforms](tutorials/transformations.md) · [Quaternions](tutorials/quaternions.md) · [Camera coordinates](tutorials/camera.md)

</div>
<div class="gem-card" markdown>

### Curves and geometry

Quadratic/cubic Bezier evaluation and adaptive sampling, Plane/Ray
representations and their supported geometric operations.

[Bezier paths](tutorials/bezier.md) · [Planes and rays](tutorials/geometry.md) · [Geometry APIs](api/plane.md)

</div>
<div class="gem-card" markdown>

### Functions and lighting

Legendre functions and real spherical harmonics: project radiance, rotate
coefficients analytically and reconstruct diffuse irradiance.

[Lighting tutorial](tutorials/lighting.md) · [Legendre](api/legendre.md) · [SH API](api/spherical-harmonics.md)

</div>

</div>

## See the mathematics

<figure markdown>

[![Row-vector transformations and composition](../examples/showcase/output/transforms.svg)](examples/gallery/transforms.md)

<figcaption markdown>
Scale, rotation and translation evaluated through gem's Matrix APIs. Dashed
shapes are the original; solid shapes show the transformed result.
[Inspect the coordinates and composition order](examples/gallery/transforms.md).
</figcaption>

</figure>

<div class="gem-comparison" markdown>

<figure markdown>

[![Diffuse sphere under the original asymmetric HDR lighting](../examples/output/sh_original.png)](examples/gallery/lighting.md)

<figcaption>Original L2 spherical-harmonic lighting.</figcaption>

</figure>
<figure markdown>

[![The same sphere after active plus-90-degree Z lighting rotation](../examples/output/sh_rotated.png)](examples/gallery/lighting.md)

<figcaption>Active +90° Z rotation; identical camera, material and display settings.</figcaption>

</figure>

</div>

These CPU reference images are generated from actual gem calculations. The
[HDR workflow](../examples/hdr_sh/README.md) connects canonical RGB coefficients
to Python and GLSL evaluation. Explore the [visual gallery](examples/index.md)
or [regenerate its assets](../examples/showcase/README.md).

## Know the conventions

Matrices use row-major storage and **row-vector mathematics**: `M * v` evaluates
the row product vM, and `A * B` applies A then B. Quaternions store **[w,x,y,z]**
and use Hamilton multiplication. Returning operations and in-place operations
have explicit ownership contracts.

Read [conventions](architecture/conventions.md) and [mathematical notation](architecture/notation.md)
before connecting another graphics API. The [numerical accuracy tutorial](tutorials/numerical.md)
explains finite precision and the limits of the supported domains.

## Choose a learning path

- **Movement and animation:** [vectors](tutorials/vectors.md) →
  [transforms](tutorials/transformations.md) → [quaternions](tutorials/quaternions.md) →
  [Bezier paths](tutorials/bezier.md).
- **Camera and picking:** [transforms](tutorials/transformations.md) →
  [camera coordinates](tutorials/camera.md) → [planes and rays](tutorials/geometry.md).
- **Lighting references:** [numerical accuracy](tutorials/numerical.md) →
  [SH lighting](tutorials/lighting.md) → [HDR environment workflow](../examples/hdr_sh/README.md).

The [tutorial index](tutorials/index.md) lists all nine tutorials, prerequisites
and independent checks. Use the [API reference](api/index.md) for exact signatures
and input/output contracts.

## Development status

gem 1.0 is in preparation with its mathematical API feature scope frozen.
Version 1.0.0 is prepared but not published; historical PyPI 0.1.12 does not contain
these audited modernizations. CPython 3.10–3.14 is tested on Linux x86_64.
Other platforms and PyPy remain unverified; Python 2.7 is unsupported.

See [compatibility](architecture/compatibility.md), [release information](development/releases.md)
and the [canonical roadmap](../ROADMAP.md) for supported,
legacy and planned functionality. gem supplies mathematical primitives for
graphics tools and reference calculations; rendering, visibility tracing and
physics systems remain separate concerns.

Created by Alex Marinescu. Source and examples use the
[BSD 2-Clause license](../LICENSE). See [contributing](development/contributing.md)
and [website verification](development/visual-overhaul.md).
