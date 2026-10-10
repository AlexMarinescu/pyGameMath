# Phase 5A — Packaging and Python compatibility

Base: `15fbce64fa87008203890142070dbf30ea802ebc` (merged PR #60).
Branch: `release/phase5a-packaging`. Prepared version: `1.0.0`.
No release, tag, upload or remote ownership change is part of this preparation.

## Summary

PEP 517/621 setuptools metadata replaces duplicated legacy setup metadata while
preserving distribution/import identity `gem`, the BSD-2-Clause license and the
sole runtime dependency `six`. `gem/_version.py` is authoritative; the initializer
only reexports `__version__`. All 16 other existing runtime files are byte-identical
to master, including compatibility shims. All 268 documented mathematical
APIs and 43 executable examples remain unchanged.

## Distribution identity and publishing authority

The public PyPI JSON record identifies `gem` 0.1.12, author Alex Marinescu,
`ale632007@gmail.com`, and the former `explosiveduck/pyGameMath` URL. Its nine
historical releases span v0.1.4–0.1.12; latest upload is June 21, 2017.
The old GitHub URL resolves to this repository's current HEAD.
The downloaded 0.1.12 archive SHA-256 is
`b236efb47897cbb37e3d07d36e2a613027e75eb96d726ad174b1ae1d270d7fd5`,
matching PyPI. All seven archived runtime Python files match repository commit
`5257291431bb45db0274dc48edf24694ecfe2e2d` after CRLF/LF normalization.
Historical distribution association is therefore established independently of
matching names or author labels alone.

**B01 — Release blocker:** current PyPI Owner/Maintainer membership and authorized
publishing credentials/trusted-publisher binding are not verified. Public author
metadata and GitHub access do not establish those permissions. The public HTML
page does not provide usable account evidence in this environment; the JSON
metadata is not an access-control list. Cloud runtime readiness reports no
configured publishing secrets or outbound publishing identities. No credential
values, private configuration files or package-upload endpoints were inspected.
An authorized maintainer must confirm PyPI membership, 2FA requirements and a
project-scoped token or approved trusted publisher before release. Keep the
existing name; there is no demonstrated naming conflict requiring a substitute.

## Build and version decisions

- `pyproject.toml` declares setuptools.build_meta, Python >=3.10, current project
  URLs, Markdown README, SPDX BSD-2-Clause and explicit LICENSE inclusion.
- Version 1.0.0 is a preparation value, not evidence of publication. The literal
  in `gem/_version.py` drives metadata and `gem.__version__`; installed
  `importlib.metadata.version("gem")` agrees. The literal is readable during an
  isolated build without importing six or runtime mathematical modules.
- setup.py remains a metadata-free compatibility entry point; setup.cfg no longer
  duplicates metadata. README.rst is retained as a pointer rather than publishing
  obsolete Python 2/PyPy and API claims.
- README links/assets use absolute repository URLs so the distribution description
  does not rely on PyPI resolving checkout-relative paths.
- Wheel naming is `gem-1.0.0-py3-none-any.whl`; source archive is
  `gem-1.0.0.tar.gz`. No entry points, numerical frameworks or compiled backends.
- Wheel: 18 Python files across `gem` and `gem.experimental`, generated metadata
  and license. Source: those files plus README, documentation, reference assets,
  tests, benchmarks/audits and verification tools. Cache/build/site/venv artifacts
  and retired transport are excluded. The final external artifact-result manifest
  is excluded from its own sdist to avoid a self-referential checksum.
- SemVer/PEP 440 policy is documented, including compatible fix/addition boundaries
  and explicit review for incompatible changes. Transitional paths have no new
  removal deadline or runtime warnings.

## Python compatibility evidence

| Interpreter | Passed | Failed/errors | Xfailed | Skipped | Installed verification |
|---|---:|---:|---:|---:|---|
| CPython 3.10.21 | 3,966 | 0 | 0 | 0 | Wheel + sdist, 268 APIs / 43 examples each |
| CPython 3.11.16 | 3,966 | 0 | 0 | 0 | Wheel + sdist, 268 APIs / 43 examples each |
| CPython 3.12.14 | 3,966 | 0 | 0 | 0 | Wheel + sdist, 268 APIs / 43 examples each |
| CPython 3.13.5 | 3,966 | 0 | 0 | 0 | Wheel + sdist, 268 APIs / 43 examples each |
| CPython 3.14.7 | 3,966 | 0 | 0 | 0 | Wheel + sdist, 268 APIs / 43 examples each |

Linux 6.18.44 x86_64, glibc 2.41; six 1.17.0 and pytest 9.1.1 throughout.
Pinned build verification: build 1.6.1, setuptools 80.9.0, wheel 0.48.0,
twine 7.0.0, pip 26.2.1. Interpreter acquisition used uv, an external development
utility; it is not a gem/build requirement. The matrix has 19,830 passing test
executions. The earlier reference total was 3,959; seven packaging-only tests
were added, and no mathematical test was removed or weakened.

Python 2.7 was investigated but no executable was available. Concrete blockers:
SH exception chaining syntax, functools.lru_cache/math.isqrt, math.isfinite and
variadic hypot. six cannot repair these. No Python 2.7 execution/compatibility
claim is made. Python 3.2–3.9 are unsupported; the >=3.10 floor aligns the core
with an actually executed modern runtime and pinned pytest matrix, rather than
claiming every earlier runtime necessarily fails every operation. CPython 3.10's
maintenance lifetime warrants a support-floor review for later releases.

**R01 — Unverified platforms:** Windows, macOS, PyPy and other architectures have
no executed matrix here. The portable wheel tag is not a platform test result.
Retain ordinary binary64/libm platform limits and native-endian raw-probe caveats.

An initial packaging-test failure on 3.13/3.14 selected `sys._base_executable`,
which lacked setuptools, despite correctly prepared virtual environments.
The test now uses the active interpreter unless GEM_BUILD_PYTHON is explicitly
set. No runtime defect was involved. Original assertions and isolated import
checks remain intact. Final complete runs above have zero failures/xfails/skips.
JUnit uses `junit_family=legacy` for the existing record_property fixture;
this avoids reporter incompatibility warnings without suppressing tests.

## Verification commands

Commands executed in `/workspace/pyGameMath` unless noted. Interpreter tools were
installed in `/tmp/phase5a-py310` through `...py314`; artifacts and clean snapshots
were outside the source checkout. Full command arrays and artifact listings are
in `phase5a-packaging-results.json` in the review checkout.

```sh
git fetch origin master
git switch -c release/phase5a-packaging origin/master
curl -fsSL https://pypi.org/pypi/gem/json -o /tmp/phase5a-pypi.json
git ls-remote https://github.com/explosiveduck/pyGameMath.git HEAD
UV_CACHE_DIR=/tmp/phase5a-uv UV_PYTHON_INSTALL_DIR=/tmp/phase5a-interpreters uv python install 3.10 3.11 3.14
```

Each tested interpreter executes, substituting 310/311/312/313/314:

```sh
/tmp/phase5a-py310/bin/python -m pytest -q -o junit_family=legacy --junitxml=/tmp/phase5a-final-py310.xml
/tmp/phase5a-py310/bin/python tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --output /tmp/phase5a-source-docs.json
```

Clean snapshots copy tracked and intended new files, excluding Git-ignored caches
and artifact trees. The first standard isolated build and pinned no-isolation
build both succeeded. Final pinned artifacts and a repeat build use:

```sh
SOURCE_DATE_EPOCH=1791593461 /tmp/phase5a-py310/bin/python -m build --no-isolation /tmp/phase5a-final-source --outdir /tmp/phase5a-final-artifacts
SOURCE_DATE_EPOCH=1791593461 /tmp/phase5a-py310/bin/python -m build --no-isolation /tmp/phase5a-repeat-source --outdir /tmp/phase5a-repeat-artifacts
/tmp/phase5a-py310/bin/python -m twine check --strict /tmp/phase5a-final-artifacts/*
/tmp/phase5a-py310/bin/python tools/verify_distribution.py --wheel /tmp/phase5a-final-artifacts/gem-1.0.0-py3-none-any.whl --sdist /tmp/phase5a-final-artifacts/gem-1.0.0.tar.gz --output /tmp/phase5a-final-inspection.json
```

Each artifact was installed in its own fresh environment for every interpreter,
with six and pinned build tools supplied locally. Example final wheel install,
executed from `/tmp` with PYTHONPATH removed:

```sh
/tmp/phase5a-310-wheel/bin/python -I -m pip install --force-reinstall --no-index --no-deps --no-build-isolation /tmp/phase5a-final-artifacts/gem-1.0.0-py3-none-any.whl
/tmp/phase5a-310-wheel/bin/python -I /workspace/pyGameMath/tools/verify_distribution.py --installed-root /tmp/phase5a-310-wheel/lib/python3.10/site-packages
/tmp/phase5a-310-wheel/bin/python -I /workspace/pyGameMath/tools/check_architecture_docs.py --examples --getting-started --api-reference --tutorials --showcase --website --package-root /tmp/phase5a-310-wheel/lib/python3.10/site-packages
```

Repeat for fresh sdist environments and each interpreter. No source-tree import
substitution is accepted. Smoke checks include exact installed runtime bytes,
version agreement, canonical/shim identity, stable Vector norms and inverse
matrices, repaired subnormal SLERP endpoints and near-pole SH values.

```sh
/tmp/phase4f-docs-env/bin/python -m mkdocs build --strict --site-dir /tmp/phase5a-site
/tmp/phase5a-py310/bin/python tools/check_site.py --site-dir /tmp/phase5a-site --output /tmp/phase5a-site-check.json
git diff --check
```

The 268 API declaration checks, seven constants/reexports and 43 executable
examples pass from source and all ten installed combinations. Historical heading
anchors and example ASTs remain checked. Website verification finds 24 notation
MathML equations and 7,400 local links/resources; strict MkDocs build passes.

## Reproducibility and safeguards

The fixed-timestamp, fixed-backend repeat compares wheel bytes and every sdist
member's content, with exact hashes/results in the machine-readable report.
In the controlled repeat, wheel bytes match exactly and all sdist member
contents match; gzip/tar timestamps make the two sdist bytes differ. No cross-platform or future-backend determinism is promised.
Expected manifest warnings for absent optional file extensions and cache-exclusion
patterns are not archive contamination or build errors.

The documentation checker freezes reviewed master, allowing only exact hashed
packaging files and two metadata-only initializer/version files. Its hard-coded
allowed set cannot exempt vector/matrix/etc. mathematics. Existing negative
scope tests remain strict. CI repeats the interpreter/artifact matrix but is not
claimed to have run remotely in this report. Required PyPI authority confirmation
and unverified platforms remain visible in the release checklist.

No release scope expansion, mathematical algorithm change or dependency removal.
