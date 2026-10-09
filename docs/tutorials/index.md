# Practical graphics mathematics with gem

Learn how to turn the existing mathematical primitives into movement, transforms,
camera coordinates, geometric queries and lighting references. Each guide explains
the calculation before checking independently known values. Complexity labels
estimate prerequisites, not runtime cost. Allow roughly 15–25 minutes per guide.

Use [current repository installation](../getting-started/installation.md), not the
historical PyPI release. All Python blocks run independently with the installed
development package, standard library and existing six dependency; no NumPy,
OpenGL context, native backend or renderer is needed. Reference verification is
CPython 3.12/Linux, not a new support matrix or Python 2.7 claim.

## Choose a learning path

| Level | Tutorial | Prerequisites / practical outcome |
|---|---|---|
| Beginner | [Vectors and movement](vectors.md) | Coordinates and basic algebra; displacement, facing and speed-limited movement |
| Beginner → intermediate | [Object transformations](transformations.md) | Vectors; transform a point group and recover local coordinates |
| Intermediate | [Camera and projection](camera.md) | Matrix composition; trace world → clip → window and back |
| Intermediate | [Quaternion orientation](quaternions.md) | Trigonometry and transforms; interpolate object/camera directions |
| Intermediate | [Bezier motion paths](bezier.md) | Vectors and a little differentiation; distinguish parameter time from distance |
| Intermediate | [Planes, rays and picking](geometry.md) | Dot products and camera tutorial; derive a local geometric query |
| Advanced | [Spherical-harmonics lighting](lighting.md) | Dot products, integration concept; project, rotate and convolve RGB lighting |
| Intermediate | [Graphics memory interoperability](interop.md) | Matrix layout, basic ctypes; inspect upload-ready memory without GPU calls |
| Beginner → advanced | [Numerical accuracy](numerical.md) | Floating-point arithmetic; choose tolerances and recognize domain limits |

For movement, follow vectors → transforms → quaternions → curves. For picking,
follow transforms → camera → geometry. For lighting references, follow vectors →
numerical accuracy → lighting → the existing [HDR/SH workflow](../../examples/hdr_sh/README.md).

## How to use the examples

Copy a complete Python block into a file or interactive environment. Blocks include
imports and assertions; displayed rounded values explain the result, while checks
use unrounded quantities. No long companion script or additional asset is required
in this phase. Local functions in examples are educational calculations, not newly
exported gem APIs. The existing HDR example remains the longer executable workflow.

[Verification and coverage](verification.md) records every checked block and the
installed-wheel/sdist runs. [Open decisions](decisions.md) links unresolved API
policies without choosing new behavior. Use the [API reference](../api/index.md)
for exact signatures, [shared conventions](../architecture/conventions.md) for
units/ownership and [ROADMAP.md](../../ROADMAP.md) for planned features.

Return to [getting started](../getting-started/README.md) or the [documentation index](../README.md).
