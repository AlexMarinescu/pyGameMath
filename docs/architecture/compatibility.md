# Python, dependency and compatibility policy

This is the current evidence and proposed release policy, not revised package
metadata. The [API inventory](api-inventory.md) describes present interfaces;
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
| CPython 3.12.14, Linux x86_64, six 1.17.0 | Phase 1–3 audits, full regression suite, clean wheel/sdist and isolated installed execution; [Phase 3E report](../../audit/PHASE3E-PERFORMANCE-GUARDS.md) | Verified reference environment; rerun on release artifacts |
| Other modern CPython versions | No executed matrix established by these reports | Candidate targets; assess and test before advertising support |
| PyPy and other implementations | No verified interpreter run | Unverified; performance and rounding may differ |
| Python 2.7 | Historical classifiers/Travis entry and retained compatibility methods, no real 2.7 execution | Deliberate legacy engineering target; currently unverified, with known blockers |
| Python 3.2–3.5 | Legacy metadata/Travis entries, no current verification | Historical claims, not current evidence of support |
| Windows, macOS, other architectures/OS versions | Historical OS classifiers; current executed evidence is Linux x86_64 | Unverified release targets, not excluded by a deliberate platform restriction |

Do not infer support for all of Python 3 from one 3.12 run. A proposed modern
matrix should cover the chosen minimum and current maintained CPython versions,
with representative Windows/macOS/Linux jobs and optionally PyPy. Select the
minimum and concrete matrix during release engineering after source/build/tooling
assessment; no minimum version is set here.

### Python 2.7 assessment

Retaining six does not make the entire present core Python-2 compatible.
`gem.spherical_harmonics` imports `functools.lru_cache`, uses `math.isqrt` and
contains `raise ... from None`; these are unavailable or invalid in Python 2.7.
Bezier sampling uses variadic `math.hypot` and `math.isfinite`, also unavailable
there. Variadic hypot and isqrt require newer Python 3 versions as well.
Division semantics and accepted scalar types need real interpreter tests, not
syntax inspection alone. None of these findings is repaired in this phase.

The pinned audit tools ([requirements-audit.txt](../../requirements-audit.txt))
use pytest 9, which requires Python >=3.10. Tests and benchmark tooling also use
modern syntax/libraries independently of core. A 2.7 validation effort would need
a suitable isolated test harness and compatible build/install tools, plus a
security/support strategy for an end-of-life interpreter. The existing setup.py
and historical Python classifiers do not establish that modern tooling installs
on 2.7. This remains a release decision with no verified 2.7 claim.

## Packaging, platform and numeric reproducibility

The distribution remains `gem`, currently declared as `v0.1.12` in
[setup.py](../../setup.py), with `gem` and transitional `gem.experimental` packages.
The repository has legacy setup.py/setup.cfg metadata, no pyproject.toml or
declared `python_requires`. Modernization, repository URLs, supported classifiers
and documentation distribution need a separate packaging phase. Markdown guides
are source-tree documentation; the current manifest explicitly includes only
selected migration guides, not the complete architecture/reference site.

The legacy [.travis.yml](../../.travis.yml) runs `launcher.py test`; that script
prints examples and does not execute pytest. It is not evidence of a current
functional CI matrix. Later CI should run real assertions and installed artifacts.

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

There is no newly declared blanket 1.0 stability guarantee. Before 1.0, review the
public surface and document its supported domains. For a future stable release,
propose semantic versioning: additive compatible interfaces in minor releases,
compatible fixes in patches, and incompatible removals/changes at a declared
major boundary. Numerical corrections still need explicit compatibility notes.
Maintainers must confirm the release policy rather than infer it from this charter.

Transitional imports retain one canonical implementation. No runtime deprecation
warnings or removal date are introduced. Retiring their paths requires a release
boundary, compatibility window and migration plan; see the
[existing migration table](../EXPERIMENTAL_MIGRATION.md). The removed E07 transport
module has no replacement. This phase changes neither shims nor licensing.
