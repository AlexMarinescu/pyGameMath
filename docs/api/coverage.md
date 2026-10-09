# API coverage and verification

## Scope and declaration coverage

[Machine-readable coverage](coverage.json) maps every audited core source
declaration to a topic page. It contains 268 entries: functions, classes, public
methods, constructors/operators and Python 3 division aliases. **None of those
declarations is deliberately omitted.** Declarations and aliases are counted
separately; coverage does not mean every possible invalid input is specified or
every method has its own example. Parameter/ownership prerequisites shared by a
group apply to each row; explicit exceptions override those group descriptions.

Seven exposed buffers/type aliases are additionally documented: common.GLfloat
and the six vector reference lists. Public instance fields are explained on their
class topic pages. Core has nine nonempty mathematical modules plus an empty
initializer. Retained reexports are covered by [legacy interfaces](legacy.md).

Excluded from the public-core count: private underscore helpers (except wrapper
special methods), incidental imported math/six/ctypes/random/etc dependencies,
source-tree launcher/examples/benchmark/doc tools, and retired E07 names. The
private historical _bezier_legacy import path remains explicitly documented.
These exclusions do not conceal a missing public declaration. `__all__` lists
in shims are export metadata, covered on the legacy page, not algorithms.

| Core module | Declarations |
|---|---|
| `gem.bezier` | 12 |
| `gem.common` | 12 |
| `gem.legendre` | 6 |
| `gem.matrix` | 57 |
| `gem.plane` | 15 |
| `gem.quaternion` | 65 |
| `gem.ray` | 7 |
| `gem.spherical_harmonics` | 18 |
| `gem.vector` | 76 |

268 declarations comprise 112 functions, 9 classes, 143 methods and 4 division aliases.

## Repeatable checks

```sh
python tools/check_architecture_docs.py --examples --getting-started --api-reference
python -m pytest -q
```

For isolated installed-wheel verification, run the checker with the fresh
installation's interpreter from outside the source tree, passing its actual
site-packages directory to `--package-root` (see
[installation verification](../getting-started/verification.md)). This mode checks
runtime signatures, constant values and reexport object identity in addition to
executing snippets. The static pass matches AST source/catalog and coverage maps,
checks Markdown links/anchors and verifies protected files against the audited base.

Examples use independent known answers/tolerances and meaningful mutation checks,
not floating-point byte identity. Coverage is a maintainability check, not a
semantic proof: source, independent regression tests and installed examples were
reviewed together. Remaining ambiguities are [documented decisions](decisions.md).

## Executed evidence

The [verification record](verification-results.json) records full-suite counts,
commands, builds, installed import paths and checker results. This phase does
not change mathematics, tests, benchmarks, runtime dependencies, packaging or
license. Wheel/sdist verification does not imply these new Markdown pages are
shipped under the unchanged distribution manifest.

Full pytest: **2,249 passed**, zero failures, expected failures or skips.
Installed wheel and independently installed sdist: **29 snippets passed**,
including 16 new API examples and 13 existing examples; 313 internal links,
268 declarations, seven buffers/type aliases and reexport identities checked.
Builds and the original Phase 4A checker passed. Missing coverage and a changed
documented signature were independently rejected by the checker.

Only CPython 3.12/Linux was executed. No new Python-version/platform support,
general invalid-input policy or release-stability claim follows from these checks.
See [index](index.md) and [compatibility](../architecture/compatibility.md).
