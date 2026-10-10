# Install current development code

Distribution: **gem**. Import namespace: **gem**. Repository: **pyGameMath**.
CPython 3.10–3.14 is tested on Linux x86_64 with six 1.17.0.
The commands below use Python 3.12 as one tested example; the declared minimum
is Python 3.10. Windows, macOS and PyPy remain unverified.

## Source installation

### Linux and macOS

Create a new environment so an older package with the same distribution/version
does not obscure the source being tested. POSIX shell:

```sh
git clone https://github.com/AlexMarinescu/pyGameMath.git
cd pyGameMath
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install .
```

### Windows PowerShell

Provided Python 3.12 and Git are installed:

```powershell
git clone https://github.com/AlexMarinescu/pyGameMath.git
Set-Location pyGameMath
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install .
```

### Windows Command Prompt

With Python 3.12 and Git installed, run these in Command Prompt:

```bat
git clone https://github.com/AlexMarinescu/pyGameMath.git
cd pyGameMath
py -3.12 -m venv .venv
.venv\Scripts\python.exe -m pip install .
```

These commands select the environment's Python directly; activation is optional.
The Windows and macOS instructions are provided for setup, not a claim that those
platforms have been tested here. Linux is the verified reference environment.

Using the environment's python executable avoids requiring activation. On Windows,
substitute `.venv\Scripts\python.exe` for `python` in subsequent commands. These
equivalents are documented, not a claim that Windows was executed in this phase.

pip installs the existing six runtime requirement; building the current source
also needs setuptools/wheel build tooling. Allow pip's normal build/dependency
resolution when using an online environment. There are no optional mathematical
or rendering dependencies needed for the quick start. Offline validation uses
locally supplied dependencies/build tools, as recorded in [verification](verification.md).

For edits to gem itself, use `python -m pip install -e .` instead of the ordinary
install. Editable imports intentionally reference the checkout; they do not prove
that a wheel works independently of the source tree.

## Verify the installed import

Run with the chosen environment's interpreter, preferably outside the checkout:

```sh
python -c "import gem; print(gem.__file__)"
python -c "from gem.vector import Vector; print((Vector(3, [1, 2, 3]) + Vector(3, [4, 5, 6])).vector)"
```

The first command identifies the installed package location; the second prints
`[5, 7, 9]`. `gem` exposes version metadata in its package initializer: import `Vector` from
`gem.vector`, not `from gem import Vector`. Follow the
[first mathematical examples](quick-start.md) once this check succeeds.

## Published package versus repository code

The [PyPI project page](https://pypi.org/project/gem/) identifies Alex Marinescu's
historical gem 0.1.12 release, uploaded June 21, 2017. The repository now prepares version 1.0.0 with the validated repairs and core
promotions. That version has not been published; the historical package is not
equivalent to current master.

The published archive was installed separately on CPython 3.12 for a basic Vector
import/addition smoke test. The canonical Bezier, Legendre and SH modules were absent;
its mathematical correctness and current feature parity are not claimed. Use source
installation for these examples rather than treating an unqualified PyPI install
as current master. Current PyPI maintainer authority remains a release blocker; no package upload
is part of release preparation.

Legacy classifiers for Python 2.7 and older Python 3 versions do not constitute
test evidence. Read the [compatibility assessment](../architecture/compatibility.md)
for known Python 2.7 blockers and the unverified platform/interpreter matrix.

## Existing build configuration

The distribution name remains `gem`; the prepared version is 1.0.0 and is not
published. Modern [pyproject.toml](../../pyproject.toml) metadata uses setuptools
and retains both the canonical modules and transitional imports. CPython
3.10–3.14 is tested on Linux x86_64. Python 2.7 and Python <3.10 are not supported;
Windows, macOS and PyPy remain unverified. `six` is retained as the only runtime
dependency and does not imply Python 2 compatibility.

The wheel contains the runtime modules, distribution metadata and BSD license.
The source distribution additionally includes documentation, tests, reference
examples and build/verification tools. See the [packaging guide](../development/packaging.md)
for exact build commands, version policy and the outstanding PyPI authority check.
