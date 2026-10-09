# Architecture documentation verification

Reference: master `dc418923c692b609a9fc611c66433e5950e0a321` (merged PR #38).
The architecture/roadmap change adds documentation and a standard-library
documentation checker only. Existing mathematical source, tests, examples,
benchmark workloads, packaging, dependencies, license, README and Wiki are unchanged.

From the repository checkout with the existing runtime/test dependencies:

```sh
python tools/check_architecture_docs.py --examples --output /tmp/architecture-docs.json
python -m pytest -q
git diff --check
```

The checker requires modern Python (AST unparse/path helpers) and the audited base
git object; it is not installed as gem runtime and makes no Python 2.7 tooling
claim. It checks local Markdown paths/heading anchors in all new navigation,
architecture and roadmap pages, matches the declaration catalog against source
and executes the four standalone Python convention examples. It also compares
protected source/test/example/benchmark/build/license files against the base.
Existing unrelated historical documentation links are not rewritten or certified.

The inventory's semantic descriptions were reviewed against implementation,
representative tests, current guides, audit findings/decisions and the historical
Wiki snapshot. Static checks verify completeness/names and executable examples;
they do not prove every prose statement or unsupported numerical domain.

Executed on CPython 3.12.14, Linux x86_64, six 1.17.0, pytest 9.1.1:
the full suite passes **2,249 tests, zero expected failures, unexpected failures
or skips**. This includes the existing clean wheel/sdist and isolated installed
checks, mathematical regressions and HDR reference tests. The documentation
checker passes 129 local links, 268 source declarations and four examples;
93 protected files are unchanged from the audited master.
These results verify this reference environment, not other advertised classifiers.

Consistency review separates implemented capabilities from provisional API sketches,
unit prerequisites from actual validation, and legacy compatibility paths from
active experimental algorithms. The [canonical roadmap](../../ROADMAP.md) is
linked rather than duplicated on the [development entry](roadmap.md). Planned
documentation sections are described without links to nonexistent pages.

Maintainer decisions are listed on the [development entry](roadmap.md#open-decisions).
No previously settled mathematical convention is reopened. Selecting support
targets, stability promises, removal windows, future interfaces or documentation/
release tooling belongs to subsequent reviewed work; none is silently implemented.
