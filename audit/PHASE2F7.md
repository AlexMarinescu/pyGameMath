# Experimental retirement and core consolidation

Base: master `62a51f9005f0cf1541f5c7996123b6b077edf847`, merged PR #31.
Branch: `fix/phase2f7-experimental-retirement`.

## Disposition and compatibility

The authoritative [Phase 2F-6 inventory](PHASE2F6.md) identifies one unfinished
implementation and seven package-marker/reexport modules. Validated Bezier,
Legendre and SH algorithms already have one canonical implementation in core;
no algorithms or numerical policies change in this phase.

`gem.experimental.sph_object` is removed, including `SPHVertex`, `SPHObject`
and `GenereateCoeffs`. Its only repository consumers were the E07 expected-failure
case and nine diagnostic tests. No core runtime or supported example imports it.
External callers are unknown; importing the retired module now fails, and old
serialized references to its containers may need application migration.
No rendering-specific scaffolding is promoted or offered as a working transport
replacement. Future scene visibility/transport remains separately designed.

The remaining experimental modules are retained as minimal reexports. They
preserve established import compatibility without algorithm duplication. Their
complete deletion conflicts with the transitional import promises in prior
phase compatibility notes. A release boundary and compatibility window must be
chosen before removing those paths. **The experimental directory is not fully
eliminated.** No runtime warnings or new APIs are added. See the public
[migration table](../docs/EXPERIMENTAL_MIGRATION.md) for each path and action.

Ordinary polynomial/SH tests now use canonical imports. Dedicated compatibility
identity tests continue to exercise all retained aliases. Core modules and
examples already use canonical imports. An AST dependency check guards against
orphaned runtime/example imports and confirms shim files contain only docstrings,
core imports and optional export lists.

## Packaging

Both core and transitional package entries remain intentional in setup.py.
MANIFEST.in includes the migration guide in source distributions. A clean build
from a temporary source copy prevents stale build/egg-info files from masking
removal. The regression builds a wheel and sdist offline, inspects archive
contents, installs the wheel with `--no-index --no-deps` into a separate target,
and runs imports and independent known answers outside the source tree under
`python -I`. Only the installed package and existing six dependency directory
are added to that interpreter. Core import origins, all alias identities and
E07 absence are verified. No mandatory dependency is added.

Development packaging tests require setuptools/wheel in the build interpreter;
`GEM_BUILD_PYTHON` selects it, defaulting to the virtualenv's base executable.
Existing legacy metadata/build warnings remain separate packaging work.

## Verification and test accounting

Environment: CPython 3.12.14, pytest 9.1.1, six 1.17.0.

```
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider --junitxml=/tmp/phase2f7.xml
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_core_packaging.py -q -p no:cacheprovider
/workspace/.venvs/pyGameMath/bin/python -m pytest tests/test_vector_common.py -k 'viewport_vector or clamp_preserves_input or empty_equality or equality_dimensions' --runxfail -q -p no:cacheprovider
```

Baseline rerun: **1780 passed, 5 xfailed**.
Final: **1778 passed, 4 xfailed, 0 unexpected failures, 0 skips**.
Packaging/architecture regressions: **7 passed**, including isolated installed
wheel and sdist checks. Deliberately retired: **10 cases** (one E07 xfail and
nine diagnostic passes). Added: seven passing cases. Thus total collected cases
change from 1785 to 1782. Archived text preserves both test bodies under
[retired-tests](retired-tests/README.md); source and implementation remain in git
history. The E07 xfail is retired explicitly, not reported as mathematically fixed.

All remaining expected failures were investigated using `--runxfail`: **4 failed,
45 deselected**. Their exact identities and observations remain unchanged:

| Test | Observation | Disposition |
| --- | --- | --- |
| `tests.test_vector_common::test_viewport_vector` | C02: getViewPort subscripts Vector, raising TypeError | Confirmed defect, separate correction |
| `tests.test_vector_common::test_equality_dimensions` | Unequal-sized vectors can compare true using smaller receiver's components | QD01 dimension policy unresolved |
| `tests.test_vector_common::test_empty_equality` | Empty equality returns None rather than proposed True | QD01 empty-domain policy unresolved |
| `tests.test_vector_common::test_clamp_preserves_input` | Returning clamp also changes caller value list | V04/QD02 ownership policy unresolved |

No remaining expected-failure marker is removed, changed or converted to a skip.
Independent existing Bezier, Legendre and SH numerical tests, ownership tests,
rotation/projection references and example tests remain in the full suite.
Machine-readable counts are in `phase2f7-test-results.json`.

## Remaining release decision

Choose the release boundary and compatibility period for the seven remaining
package-marker/reexport files. Until that decision, keep the tiny alias layer
and its wheel/identity tests. At final removal, update package declarations,
public migration status and import consumers together, clean build artifacts,
and verify absence from source distributions and installed wheels. No Phase 3
optimization, release publication or merge is included.
