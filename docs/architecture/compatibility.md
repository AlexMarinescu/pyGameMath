# Python, dependency and compatibility policy

This is the current release-preparation evidence and compatibility policy. The [API inventory](api-inventory.md) describes present interfaces;
the [roadmap](../../ROADMAP.md) separates completed work from release preparation.

## Runtime boundary

Mathematical implementations stay pure Python. Standard-library `math`, `ctypes`,
`random`, `struct` and other utilities are allowed. NumPy, Cython, compiled
mathematical extensions, native SIMD backends and mandatory numerical frameworks
are excluded. The existing runtime dependency is `six`; its compatibility helpers
remain in use. Reducing runtime dependencies is a long-term objective, not
permission to remove six before evaluating interpreter compatibility.

`ctypes` is standard-library interoperability, not an accelerated mathematical
backend. `Matrix.c_matrix` is a row-major float32 snapshot of Python storage.
Constructor and supported in-place matrix operations refresh it; direct external
list edits do not. Do not change its type, order or ownership as an optimization.
Foreign rendering consumers must explicitly reconcile their upload layout and
transpose conventions; no OpenGL driver support claim follows from a ctypes test.

## Verified evidence and unverified targets

| Interpreter/platform | Evidence | Policy status |
|---|---|---|
| CPython 3.10.21, 3.11.16, 3.12.14, 3.13.5, 3.14.7; Linux x86_64; six 1.17.0 | Phase 5A full suite and installed wheel/sdist checks | Tested modern release matrix; keep testing maintained patch versions |
| PyPy and other implementations | No executed interpreter run | Unverified |
| Python 2.7 | Source contains unsupported syntax and standard-library calls; no interpreter available | Unsupported; classifiers removed, no compatibility rewrite |
| Python 3.2–3.9 | Historical metadata is not evidence; minimum is now 3.10 | Unsupported |
| Windows, macOS and other architectures | No executed platform matrix | Unverified; pure-Python wheel is not proof of platform testing |

### Python 2.7 assessment

Retaining six does not make the core Python-2 compatible. SH contains
`raise ... from None`, imports `functools.lru_cache`, and calls `math.isqrt`.
Bezier uses variadic `math.hypot` and `math.isfinite`. These have concrete
Python 2.7 syntax or library incompatibilities. No Python 2.7 interpreter was
available for execution, and no compatibility pass is claimed.

The pinned audit tooling uses pytest 9 (Python >=3.10). The chosen 3.10 minimum
provides a single executed runtime/test matrix without rewriting mathematics to
serve end-of-life interpreters. Python 3.8/3.9 may support some core functions but
are outside the declared and tested range; that is a support-policy boundary,
not a claim that every import necessarily fails there.

## Packaging, platform and numeric reproducibility

The distribution and import namespace remain `gem`. Version 1.0.0 has one source,
`gem._version.__version__`, and is exposed as `gem.__version__`. Metadata is in
[pyproject.toml](../../pyproject.toml); setup.py is a compatibility entry point.
The only runtime dependency remains `six`. The pure-Python wheel includes all
core modules and existing transitional `gem.experimental` reexports. The sdist
includes documentation and verification sources. See the
[packaging policy](../development/packaging.md) for contents and release blockers.

The historical [.travis.yml](../../.travis.yml) is not an active support matrix:
its `launcher.py test` prints examples and does not run pytest. Current packaging
checks use real pytest, built artifacts and isolated imports; the new GitHub
workflow is a repeatable check, not evidence that a remote CI run has completed.

Python float calculations ordinarily use binary64; caller numeric types are not
uniformly coerced, and no arbitrary-precision public contract is established.
ctypes exports use binary32 and may round or overflow independently of Python
results. Raw angular probes are native-endian float32 RGB: convert endianness
explicitly when exchanging files across platforms.

Use independent references and operation-specific tolerances; do not impose a
single epsilon on every domain. Exact Vector equality stays exact. Matrix inverse
singularity detection has no arbitrary near-singular threshold. Stable finite
norms/inverses do not imply stable determinants, dot products or every quaternion
operation at extreme scales. Numerical, integration and truncation limits are
described in [conventions](conventions.md).

libm/interpreter rounding and zlib output can differ across systems. HDR image
hashes are a same-environment reproducibility check; coefficient, linear-pixel
and mathematical-reference tests remain the portable oracle. Performance captures
from different environments are not controlled comparisons.

## API preservation and versioning

Preserve current names, signatures, return types, caller ownership, receiver
mutation and mathematical conventions. Float-only Matrix/Quaternion wrapper
division, mixed angle units, misspelled methods and legacy forward axes remain
intentional compatibility obligations until explicitly reviewed. Returning methods
usually allocate fresh storage; constructors and explicit accessors are exceptions.

Version 1.0.0 is prepared for review, not published. The release versioning
policy uses semantic versioning: additive compatible interfaces in minor releases,
compatible fixes in patches, and incompatible removals/changes at a declared
major boundary. Numerical corrections still need explicit compatibility notes.
Unsupported domains do not acquire new guarantees from the version number.

Transitional imports retain one canonical implementation. No runtime deprecation
warnings or removal date are introduced. Retiring their paths requires a release
boundary, compatibility window and migration plan; see the
[existing migration table](../EXPERIMENTAL_MIGRATION.md). The removed E07 transport
module has no replacement. This phase changes neither shims nor licensing.
