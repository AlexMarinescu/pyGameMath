# gem 1.0.0 release notes (prepared)

Pure-Python graphics mathematics: vectors, matrices, quaternions, planes, rays,
Bezier curves, Legendre polynomials and spherical harmonics. The gem 1.0
mathematical API is frozen, with 268 documented declarations and 43 executable
examples. Matrices use row-major storage and row-vector transformations;
quaternions use [w,x,y,z] and Hamilton multiplication.

The release strategy is GitHub-first. This document prepares release notes;
it does not establish that a release has been published. PyPI 1.0.0 publication
is deferred until publishing authority is verified.

CPython 3.10–3.14 are the supported interpreter targets. Completed CI runs must
establish platform verification before release; workflow configuration alone
is not evidence of Windows/macOS support. Python 2.7 is unsupported.

Download the approved wheel or sdist, verify SHA-256 checksums, then install
with `python -m pip install ./gem-1.0.0-py3-none-any.whl` or the corresponding
sdist path. Do not use `pip install gem==1.0.0` as GitHub installation guidance.
The only mandatory runtime dependency is six. No OpenGL context or compiled
numerical backend is required. Optional example/documentation dependencies
are separate.

Historical experimental import shims remain transitional; canonical modules
are gem.bezier, gem.legendre and gem.spherical_harmonics. Unfinished shadow
transport is retired. Numerical guarantees remain subject to binary64 limits
and documented domains; SH lighting is a low-order approximation, not a
visibility or rendering engine. See the mathematical conventions and
compatibility documentation for detailed contracts.

Preserve the approved artifacts unchanged for eventual PyPI upload. Do not
rebuild or overwrite official release assets silently.
