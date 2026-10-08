# Phase 2B matrix division and inversion consistency

Branch: `fix/phase2b-matrix-division`, based on latest master `c05e9cc973eee1614c9d0537ca75a99e0aad52cf`, the merge of PR #11. Scope is M02, M03, M07 and their inverse4 dependency from [PHASE2-DECISIONS.md](PHASE2-DECISIONS.md). Phase 2A and historical audit evidence are preserved.

## Implementation and compatibility

`matrix_div` now divides elements without exchanging row and column indices. `inverse4` explicitly transposes its cofactor matrix to obtain the adjugate before division. Python 3 division special methods expose the existing legacy implementations; both legacy names remain. In-place division refreshes `c_matrix` and returns self.

The wrapper's existing float-only scalar domain is preserved; adding integer or broader numeric protocol support remains QD12. Storage, row-vector application, matrix multiplication, receiver ownership, singular/zero exceptions, and pure-Python dependencies are preserved. Numerical division results change for non-symmetric matrices, and in-place ctypes snapshots now show the divided values. Migration examples and detailed implications appear in [COMPATIBILITY.md](COMPATIBILITY.md#phase-2b-m02-m03-m07).

## Validation

CPython 3.12.14, pytest 9.1.1, six 1.17.0. The unchanged merged baseline reproduced **312 passed / 69 xfailed**, exit 0. Before corrections, a focused `--runxfail` run of all four original M02/M03/M07 cases produced **4 failed**, exit 1.

Full suite after corrections:

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/phase2b-results.xml
```

**421 cases: 356 passed, 65 strict xfails, no failures, XPASS, or skips; exit 0.** Four original xfails are now ordinary passes (M02: one, M03: two, M07: one), and 40 new cases pass. Remaining xfails comprise **57 other confirmed-defect cases and eight contract questions**. Comparing JUnit case identities against the baseline verifies that exactly the four targeted cases were removed from the expected-failure set; every unrelated marked case remains.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy \
  --junitxml=/tmp/phase2b-runxfail.xml
```

**356 passed / 65 failed**, expected exit 1. The failures match the remaining expected-failure set exactly. The library still has unresolved defects outside this phase. `git diff --check` passes. [Machine-readable results](phase2b-test-results.json) record counts by module and finding.

New regressions cover:

- Non-symmetric 2x2, 3x3, and 4x4 division by positive, negative, and fractional float scalars: helper, `/`, legacy `__div__`, `/=`, and legacy `__idiv__`; element positions, result types, independent rows, input preservation, receiver identity, and export synchronization.
- Known-answer non-symmetric inverses at all three sizes, including a 4x4 triangular inverse with off-diagonal entries `-2,6,-24` that distinguish the adjugate from its transpose. Both multiplication orders are checked for original and divided matrices; an independent Fraction Gauss-Jordan oracle and in-place export checks supplement known answers.
- Zero-division state preservation at all three sizes for all four division special methods; unsupported operands preserve NotImplemented delegation and receiver state.

Existing 60 seeded 2x2/3x3/4x4 inverse-oracle cases and both-sided identity checks remain ordinary passes, as do existing row-vector storage/composition and in-place method tests. Existing extreme-scale inverse4 expected failures remain. Only Python 3.12 was executed; no interpreter support policy was changed.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/matrix.py` | Correct elementwise division, explicit inverse4 adjugate, Python 3 methods, export refresh |
| `tests/test_matrix.py` | Convert M02/M03 regressions to ordinary passing tests |
| `tests/test_api_edges.py` | Convert M07 ctypes regression to an ordinary passing test |
| `tests/test_wiki_contracts.py` | Convert the documented Python 3 division example to an ordinary passing test |
| `tests/test_matrix_division.py` | Forty division, inverse, ownership, export, and error regressions |
| `audit/COMPATIBILITY.md` | Behavioral implications, retained operand domain, and migration guidance |
| `audit/PHASE2B.md` | Scope, verification, and complete file summary |
| `audit/phase2b-test-results.json` | Machine-readable verification summary |

No unrelated implementation changes, further Phase 2 work, merge into master, or PyPI publication are included.
