# pyGameMath · gem

**Math tools for Python games, graphics and simulations.**

gem helps you calculate movement, rotate objects, transform 3D coordinates,
create smooth curves and approximate lighting—all using regular Python.
No NumPy or compiled extensions are required. Its only runtime dependency is
[six](https://pypi.org/project/six/).

## See what it does

[![Sphere lit from the right](examples/output/sh_original.png)](docs/examples/gallery/lighting.md)
[![Same sphere with lighting rotated upward](examples/output/sh_rotated.png)](docs/examples/gallery/lighting.md)

**Rotate the lighting, keep the object still.** These images show the same sphere
before and after a 90-degree lighting rotation. They are calculated on the CPU;
you don't need a graphics card or an OpenGL window to run the example.
[Spherical harmonics](docs/tutorials/lighting.md) compress light arriving from many
directions into a small set of numbers, making soft lighting easier to calculate.
[Try the lighting example](examples/hdr_sh/README.md).

<a href="docs/examples/gallery/vectors.md"><img src="examples/showcase/output/vectors.svg" width="400" alt="Vector arrows showing directions and addition"></a>
<a href="docs/examples/gallery/transforms.md"><img src="examples/showcase/output/transforms.svg" width="400" alt="A rectangle moved, rotated and resized in different orders"></a>

**Directions and movement** — [add vectors](docs/tutorials/vectors.md) to combine
movement. **Object placement** — [use matrices](docs/tutorials/transformations.md)
to move, rotate and resize a shape. The order of those steps matters.

<a href="docs/examples/gallery/quaternions.md"><img src="examples/showcase/output/quaternions.svg" width="400" alt="Coordinate axes turning smoothly between two orientations"></a>
<a href="docs/examples/gallery/bezier.md"><img src="examples/showcase/output/bezier.svg" width="400" alt="Smooth Bezier curves with control points and sample markers"></a>

**Smooth turns** — [quaternions](docs/tutorials/quaternions.md) describe rotations
and help blend between them. **Smooth paths** — [Bezier curves](docs/tutorials/bezier.md)
let you shape a curve with a few control points. Select a diagram for a larger
view, its explanation and runnable source.

Every image above comes from the committed, reproducible examples.
[Explore the gallery](docs/examples/index.md) or
[regenerate the diagrams](examples/showcase/README.md).

## Install the development version

For the current features, install from this repository. **The older package on
PyPI does not include these updates.** In a terminal with Git and Python 3.12:

```sh
git clone https://github.com/AlexMarinescu/pyGameMath.git
cd pyGameMath
python3.12 -m venv .venv
.venv/bin/python -m pip install .
```

These commands are for Linux and macOS. On Windows, and for help setting up a
separate Python environment, follow the
[installation guide](docs/getting-started/installation.md). It covers PowerShell,
Command Prompt, development edits and package details.

## Try a simple movement calculation

A vector is a list of numbers representing a position, direction or movement.
Here, three numbers give the X, Y and Z coordinates. Add a movement to a position
to find where it ends up:

```python
from gem.vector import Vector

position = Vector(3, [1, 2, 3])
movement = Vector(3, [4, 0, -1])
new_position = position + movement
print(new_position.vector)  # [5, 2, 2]
```

The original position stays unchanged. To run this with the environment above,
save it as `move.py` and use `.venv/bin/python move.py`.
The [quick start](docs/getting-started/quick-start.md) also covers matrices,
rotations, curves and lighting.

## Tools available today

- **[Vectors](docs/api/vector.md):** calculate directions, distances and movement.
- **[Matrices](docs/api/matrix.md):** move, rotate and resize objects in 2D or 3D;
  convert coordinates for a camera or screen.
- **[Quaternions](docs/api/quaternion.md):** rotate objects and blend smoothly
  between orientations.
- **[Bezier curves](docs/api/bezier.md):** build smooth paths, evaluate points
  along them and sample curved sections more closely.
- **[Planes](docs/api/plane.md) and [rays](docs/api/ray.md):** describe flat
  surfaces and directed lines, and move or rotate them. Ray intersection queries
  are not implemented.
- **[Legendre functions](docs/api/legendre.md) and
  [spherical harmonics](docs/api/spherical-harmonics.md):** provide the building
  blocks for approximating light from the surrounding environment, including
  rotating that lighting and calculating its effect on a diffuse surface.

These are math tools, not a game engine or renderer. You can use them on their own
or connect them to your graphics application. The examples explain the coordinate
and rotation rules before showing how to combine operations.

## Find your next step

- [Getting started](docs/getting-started/README.md) — installation and first examples.
- [Tutorials](docs/tutorials/index.md) — learn through practical calculations.
- [API reference](docs/api/index.md) — functions, arguments and return values.
- [Visual examples](docs/examples/index.md) — diagrams and reproducible lighting.
- [Documentation home](docs/index.md) — browse the complete documentation.
- [Build the documentation website locally](docs/development/website.md) — the
  website is prepared, but no public deployment is claimed.
- [Contributing](docs/development/contributing.md) — changes, tests and review.

## What's planned?

The library is being prepared for a **1.0 release; it is not released yet**.
Later work is proposed for finding relationships between shapes, generating
repeatable noise for procedural content, measuring distance to surfaces (signed
distance fields), working with 3D grids (voxels), and calculating light through
volumes such as fog. Future C and C++ versions are proposed as separate projects.

These capabilities are **planned, not available today**, with no promised dates.
The [development roadmap](ROADMAP.md) gives the full direction and priorities.

## Compatibility and release notes

The verified development environment is **CPython 3.12 on Linux**. Other Python
versions and platforms need testing; Python 2.7 support has not been verified.
The repository still declares version 0.1.12, so matching the old
[PyPI version](https://pypi.org/project/gem/) does not mean you have the same code.
A 1.0 release or a blanket promise of API stability is not claimed.

Some constructors keep references to lists you provide. Read the
[ownership and compatibility guide](docs/architecture/compatibility.md) before
sharing mutable data. The [conventions](docs/architecture/conventions.md) explain
coordinate and angle rules. Use `gem.bezier`, `gem.legendre` and
`gem.spherical_harmonics` for current code; older experimental imports remain
available through [compatibility re-exports](docs/EXPERIMENTAL_MIGRATION.md).

## License and author

Created by **Alex Marinescu**, originally for learning graphics mathematics and
personal OpenGL projects. Licensed under the [BSD 2-Clause license](LICENSE),
copyright 2015–2026 Alex Marinescu. Historical attribution remains intact.
