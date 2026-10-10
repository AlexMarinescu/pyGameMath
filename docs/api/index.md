# gem API reference

This reference documents existing public core declarations and retained compatibility
imports at merged master `3e714fe949b7a6b7724d5c0da3395ee92483265f` (PR #59).
It covers the [268-declaration inventory](../architecture/api-inventory.md), exposed
reference buffers/type aliases and public instance fields. Tables show **exact
source signatures**, including `self`; omit `self` when calling a bound method.
Alias rows identify actual targets rather than inventing annotated signatures.

| Topic | Reference |
|---|---|
| Components, arithmetic, dot/cross, norms, comparison, ownership | [Vector](vector.md) |
| Matrix2/3/4 storage, products, division, determinant, transpose, inverse, ctypes | [Matrix](matrix.md) |
| Rotation, pivot, scale, translation, shear and local homogeneous promotion | [Transformations](transformations.md) |
| Cameras, perspective/orthographic matrices, project/unproject | [Projection](projection.md) |
| Hamilton algebra, rotation, conversions, interpolation, powers/logarithms | [Quaternion](quaternion.md) |
| Evaluation, path builders, adaptive subdivision and tolerance | [Bezier](bezier.md) |
| Unnormalized associated functions and mutable recurrence scratch | [Legendre](legendre.md) |
| Real basis, projection, sampling, analytical rotation and radiometry | [Spherical harmonics](spherical-harmonics.md) |
| Coefficients, incidence, orientation, normalization and polygon fits | [Plane](plane.md) |
| Constructor ownership, stored state, duplication and rigid transforms | [Ray](ray.md) |
| Viewport, angle conversions, scalar/raw utilities and ctypes | [Common utilities](common.md) |
| Retained import paths, legacy/private/retired distinctions | [Compatibility interfaces](legacy.md) |
| Unresolved behavior requiring separate maintainer review | [Decisions](decisions.md) |
| Declaration/constant coverage, executable examples and evidence | [Coverage and verification](coverage.md) |

Start with [installation and quick start](../getting-started/README.md) if new to
gem. Import from modules, not the empty `gem` initializer. Matrix2/3/4 and
Vector2/3/4 denote dimensions of `Matrix`/`Vector`, not separate classes.

Shared [conventions](../architecture/conventions.md),
[architecture charter](../architecture/philosophy.md),
[compatibility policy](../architecture/compatibility.md) and root
[ROADMAP.md](../../ROADMAP.md) remain canonical. Topic pages explain their application,
not a second policy. Future modules remain planned rather than importable.

Core algorithms remain pure Python with six and standard-library integration.
No optional renderer/native framework is needed for examples; verified interpreter
is CPython 3.12/Linux, not a newly declared support matrix or Python 2.7 guarantee.
The reference is source-tree Markdown; packaging/framework changes are separate.
