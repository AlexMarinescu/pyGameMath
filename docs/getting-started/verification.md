# Example verification

The README and quick-start Python snippets are executable checks, not pseudocode.
The documentation checker checks local links, the architecture API inventory and
unchanged mathematical, test, benchmark, packaging and license files.

## Repeat the checks

From a checkout with the audit requirements installed:

```sh
python -m pytest -q
python tools/check_architecture_docs.py --examples --getting-started
```

To verify snippets independently of source-tree imports, build in a clean copy,
install its wheel into a fresh virtual environment with six, then run from outside
that checkout. Pass the installed environment's site-packages directory explicitly:

```sh
python -I /path/to/pyGameMath/tools/check_architecture_docs.py --examples --getting-started --package-root /path/to/venv/lib/python3.12/site-packages
```

The path placeholders must be replaced with the actual checkout and environment.
This command checks both Phase 4A examples and the new getting-started examples.
The ordinary checker defaults to source-tree imports; the explicit installed mode
records the actual gem import path in its JSON output.

## Reference run

Verification used CPython 3.12.14 on Linux x86_64, six 1.17.0, setuptools 84.0.0,
wheel 0.48.0 and pytest 9.1.1. Build tooling is separate from gem's runtime dependency.
The source and installed-wheel runs, full regression results and commands are
recorded in [the verification record](verification-results.json).

The historical PyPI gem 0.1.12 source archive was also installed in a separate
fresh environment. Its basic Vector addition produced `[5, 7, 9]`; the canonical
`gem.bezier`, `gem.legendre` and `gem.spherical_harmonics` modules were absent.
This limited smoke test does not validate historical mathematical correctness or
feature parity with the repository. The archive SHA-256 was
`b236efb47897cbb37e3d07d36e2a613027e75eb96d726ad174b1ae1d270d7fd5`.

Windows commands are documented equivalents, not executed checks. Python 2.7 and
other interpreter/platform combinations remain unverified. No release/support
matrix or packaging changes are implied by these examples.
