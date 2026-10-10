# gem 1.0.0 release notes (prepared)

Pure-Python graphics mathematics: vectors, matrices, quaternions, planes, rays,
Bezier curves, Legendre polynomials and spherical harmonics. The gem 1.0
mathematical API is frozen, with 268 documented declarations and 43 executable
examples. Matrices use row-major storage and row-vector transformations;
quaternions use [w,x,y,z] and Hamilton multiplication.

The release strategy is GitHub-first. This document prepares release notes;
it does not establish that a release has been published. PyPI 1.0.0 publication
is deferred until publishing authority is verified.

CPython 3.10–3.14 passed the Phase 5C regression and installation matrix on
GitHub-hosted Linux x86-64, Windows AMD64 and macOS ARM64. This evidence covers
the tested interpreter builds and runners, not every operating-system release
or architecture. Final merged-master validation and owner approval remain
required before publication. Python 2.7 is unsupported.

Download the approved wheel or sdist, verify SHA-256 checksums, then install
with `python -m pip install ./gem-1.0.0-py3-none-any.whl` or the corresponding
sdist path. Do not use `pip install gem==1.0.0` as GitHub installation guidance.
The only mandatory runtime dependency is six. No OpenGL context or compiled
numerical backend is required. Optional example/documentation dependencies
are separate.

Historical experimental import shims remain transitional; canonical modules
are gem.bezier, gem.legendre and gem.spherical_harmonics. Unfinished shadow
transport is retired. The showcase retains byte-exact SVGs and checks eight audited full-precision
quaternion JSON measurements using independently bounded platform rounding.
No runtime mathematics or committed reference image is changed. Numerical
guarantees remain subject to binary64 limits
and documented domains; SH lighting is a low-order approximation, not a
visibility or rendering engine. See the mathematical conventions and
compatibility documentation for detailed contracts.

Preserve the approved artifacts unchanged for eventual PyPI upload. Do not
rebuild or overwrite official release assets silently.


GitHub Release title: **gem 1.0.0**. After publication, obtain assets from
<https://github.com/AlexMarinescu/pyGameMath/releases/tag/v1.0.0>. That URL is the
planned release location; these notes do not claim it currently exists.
Verify the downloaded `SHA256SUMS` before installation:

```sh
sha256sum -c SHA256SUMS
python -m pip install ./gem-1.0.0-py3-none-any.whl
```

On Windows, use `Get-FileHash -Algorithm SHA256` for the wheel and sdist and
compare both values with the approved checksum listing. Checksums alone do
not establish authenticity; obtain them from the reviewed release.
