# Phase 2A core mathematical correctness

Branch: `fix/phase2a-core-math`, based on latest master `75c3fde1541007f9803c903674525a3e0bb49384` (merged PRs #9 and #10). Read PHASE1, PHASE1B-WIKI, CONVENTIONS, COMPATIBILITY, and PHASE2-DECISIONS before implementation. Scope stops at V01, M01, Q01.

## Changes and compatibility

V01 compares every component of equal-size, nonempty vectors with exact numeric comparisons. M01 fills both rows with the correct adjugate divided by determinant. Q01 negates the three imaginary components before division by squared magnitude, preserving `[w,x,y,z]`. Full compatibility implications and concrete before/after examples are in [COMPATIBILITY.md](COMPATIBILITY.md#phase-2a-v01-m01-q01).

No public signatures, namespace, storage or multiplication conventions, dependencies, unrelated algorithms, packaging, or version policy change. Empty-vector behavior and unresolved dimension contracts remain deferred. Singular matrix and zero-quaternion inverse exception behavior is retained. No merge or PyPI publication is included.

## Verification

CPython 3.12.14, pytest 9.1.1, six 1.17.0. The unchanged merged baseline reproduced **254 passed / 94 strict xfails**, exit 0. Focused `--runxfail` reproductions for the original equality cases, inverse2 known answer, and quaternion inverse cases produced **5 failed / 1 passed**, exit 1, before corrections.

After corrections:

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/phase2a-results.xml
```

**381 cases: 312 passed, 69 strict xfails, no failures, no XPASS, no skips; exit 0.** The remaining xfails comprise 61 confirmed-defect cases and eight contract questions. All 25 targeted xfails become ordinary passing tests: V01 two cases, M01 21 cases, Q01 two cases. An additional 33 cases pass. Existing SyntaxWarnings in unrelated algorithms remain; seven warnings were emitted when changed modules were recompiled.

The full `--runxfail` run (`--tb=no`, same cache/JUnit settings, `/tmp/phase2a-runxfail.xml`) produced **312 passed / 69 failed**, expected exit 1. Its failures correspond exactly to the remaining marked defects/questions. This is evidence that Phase 2A did not resolve the entire library audit.

Regression coverage:

- V01 checks equality and inequality, identical vectors, differences at each component in dimensions 1,2,3,4,8, both operand orders, a 1e-12 difference to reject approximate comparison, and unsupported-operand delegation.
- M01 checks four known answers (positive/negative determinant, diagonal, zero diagonal), helper and both Matrix inverse methods, both multiplication orders against identity, input preservation, and ctypes export. Existing 20 seeded 2x2 cases compare against an independent Fraction Gauss-Jordan oracle and both identities. Existing 3x3/4x4 tests remain ordinary passes.
- Q01 checks seven known answers (unit, nonunit, real, each imaginary axis, mixed signs), helper and Quaternion method, both Hamilton multiplication orders against identity, fresh result and input preservation. Twenty seeded ordinary nonunit quaternions supplement known answers.

`git diff --check` passes. Machine-readable outcomes are in [phase2a-test-results.json](phase2a-test-results.json). Only Python 3.12 was executed; supported interpreter policy is unchanged.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/vector.py` | V01 equality/inequality loop correction |
| `gem/matrix.py` | M01 2x2 inverse row/sign correction |
| `gem/quaternion.py` | Q01 conjugation signs in inverse |
| `tests/test_vector_common.py` | Remove V01 defect marker; exact all-component regressions |
| `tests/test_matrix.py` | Remove M01 markers; known answers, both identities and wrapper checks |
| `tests/test_quaternion.py` | Remove Q01 markers; known answers and both-sided nonunit identities |
| `audit/COMPATIBILITY.md` | Phase 2A behavioral implications and deferred contracts |
| `audit/PHASE2A.md` | Scope, verification, and complete file summary |
| `audit/phase2a-test-results.json` | Machine-readable final test outcomes |

Historical audit results and wiki snapshots are preserved. No further Phase 2 fixes are started.
