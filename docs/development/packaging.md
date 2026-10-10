# Packaging and release preparation

## Identity and release status

The distribution is **gem**, the import namespace is **gem**, and the repository
is **pyGameMath**. Version **1.0.0** is prepared for review and has not been
published. No release tag or publishing workflow is introduced.

PyPI's historical 0.1.12 metadata matches Alex Marinescu, the historical email
and the former `explosiveduck/pyGameMath` repository URL. The downloaded archive
matches PyPI's SHA-256 and the repository's historical modules. This establishes
association, not current account authority. **Before release, an authorized
maintainer must verify current PyPI Owner/Maintainer membership, two-factor
requirements and a project-scoped token or approved trusted-publisher binding.**
No secrets were read, credentials tested by upload, ownership changed or
alternative distribution name selected. Public metadata cannot establish access
rights. See the [Phase 5A report](../../audit/PHASE5A-PACKAGING.md).

## Version source and policy

`gem/_version.py` is the only version declaration. `gem.__version__` reexports it;
setuptools reads the literal without importing the mathematical modules.
`importlib.metadata.version("gem")` agrees for an installed distribution.
The initializer does not reexport mathematical classes.

Use PEP 440 versions and semantic versioning: patches for compatible fixes,
minor versions for compatible additions, major versions for incompatible public
changes. Numerical fixes require compatibility notes even when signatures stay
the same. Prerelease candidates can use `1.0.0rcN` by editing the same source;
this preparation uses `1.0.0` and does not itself authorize a release.
Artifacts are `gem-1.0.0-py3-none-any.whl` and `gem-1.0.0.tar.gz`.
Existing compatibility imports have no new removal date or runtime warning.

## Tested support matrix

CPython 3.10.21, 3.11.16, 3.12.14, 3.13.5 and 3.14.7 are exercised on Linux
x86_64. The declared minimum is 3.10. Run current maintained patch versions in
CI and before release; this report is evidence for the exact versions listed,
not every patch/platform combination. Windows, macOS, PyPy and other architectures
are unverified. Python 2.7 and 3.2–3.9 are outside the supported range.

Python 2.7 has `raise ... from None` syntax and `lru_cache`, `isfinite`, `isqrt`
and variadic `hypot` incompatibilities in current core. No real 2.7 interpreter
was available; no compatibility pass is claimed. `six` remains the sole runtime
dependency for existing helpers, not a declaration of Python 2 support.

## Build from a clean checkout

```sh
python3.12 -m venv .venv-build
.venv-build/bin/python -m pip install -r requirements-build.txt
.venv-build/bin/python -m build
.venv-build/bin/python -m twine check --strict dist/*
.venv-build/bin/python tools/verify_distribution.py   --wheel dist/gem-1.0.0-py3-none-any.whl --sdist dist/gem-1.0.0.tar.gz
```

The standard build command uses PEP 517 isolation and builds the wheel from the
sdist. For calibrated/offline reproduction, preinstall the pinned build tools,
set `SOURCE_DATE_EPOCH` to the base commit timestamp, and use `python -m build
--no-isolation`. Pinning verification tools records this run; the build backend
requirement in pyproject.toml permits compatible newer setuptools releases.
`setup.py` remains a minimal compatibility entry point; direct invocation is
legacy, not the recommended release workflow.

The wheel contains only the two runtime packages, generated metadata and license;
there are no command-line entry points. The sdist includes core, license, README,
documentation, examples/assets, mathematical tests, benchmark/audit records and
verification tooling. These are deliberate source contents, not runtime dependencies.
Caches, build trees, virtual environments and retired shadow transport are excluded.
The complete MkDocs site still needs its optional documentation requirements.
Git-history-based scope checks require a checkout, not just the sdist.

## Installed verification

Create separate environments for the wheel and sdist. Install each with pip from
outside the checkout, then verify the metadata and source hashes:

```sh
python3.12 -m venv /tmp/gem-wheel
/tmp/gem-wheel/bin/python -m pip install dist/gem-1.0.0-py3-none-any.whl
/tmp/gem-wheel/bin/python -I tools/verify_distribution.py   --installed-root /tmp/gem-wheel/lib/python3.12/site-packages
```

Repeat with a fresh environment and the tarball. For an offline test, preinstall
`six` and the build backend there, then use `--no-index --no-deps
--no-build-isolation`. The report includes the exact commands actually run.
The API/example checker executes all 268 declarations and 43 documentation
examples against an explicit installed root; it also preserves original example
ASTs and heading anchors. Packaging scope changes are accepted only by reviewed
content fingerprints; mathematical files remain frozen.

## Release checklist

- Confirm the current PyPI account authority and publishing configuration.
- Review this PR, version, support matrix, license and source contents.
- Run the full suite on maintained interpreter patches and supported platforms.
- Build from the reviewed commit; validate metadata, imports and runtime hashes.
- Execute all documentation examples against wheel and sdist installs.
- Review artifact hashes and compatibility notes; retain provenance.
- Require a separate human-controlled publishing/tagging step.

Outstanding risks include unverified non-Linux platforms/PyPy, historical import
transition obligations and binary64/libm platform rounding. Matching artifacts
under fixed tools/timestamps is a reproducibility check, not a guarantee across
backend releases or operating systems.
