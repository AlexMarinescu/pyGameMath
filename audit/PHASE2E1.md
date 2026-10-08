# Quaternion/matrix conversion correctness

Branch `fix/phase2e1-quaternion-matrix` starts from master `b5a01e9dfaaa469365f1894538db614a129d33ea`, the merge of PR #17. This change addresses Q02/Q07.

## Changes

Matrix-to-quaternion conversion updates the largest squared-component candidate through all four branches and reconstructs the quaternion from that value. Quaternion-to-matrix conversion refreshes the completed matrix's float32 ctypes snapshot. Quaternion order [w,x,y,z], row-vector rotation signs, unit-rotation prerequisites, signatures, return types, and input ownership are preserved.

Historical source `4253839` contains the same selection defects; the quaternion wiki is a placeholder. Existing source and audit tests establish rotation orientation and sign equivalence. [CONVENTIONS.md](CONVENTIONS.md#quaternionmatrix-conversions) documents the conversion domain; [COMPATIBILITY.md](COMPATIBILITY.md#quaternionmatrix-conversions) describes corrected values and exports. No quaternion normalization, invalid-matrix validation, forward-axis change, or unrelated algorithm correction is introduced.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Unchanged master: **800 passed, 48 strict xfails**. The four Q02 half-turn cases reproduce ZeroDivisionError/ValueError under `--runxfail`, and Q07 reproduces the stale identity export. The 67 new cases produced **62 failed, five passed** before correction.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/quat-matrix-results.xml
```

Full suite: **915 cases, 872 passed, 0 failed, 43 strict xfails**, exit 0, with no XPASS or skips. Converts four Q02 cases and one Q07 case and adds 67 passing regressions.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/quat-matrix-runxfail.xml
```

This run produced **872 passed, 43 failed**, expected exit 1. JUnit case comparison verifies exact preservation of every unrelated expected failure, with exposed failures matching the remaining marker set. `git diff --check` passes. [Machine-readable results](phase2e1-test-results.json) record module/finding counts and all remaining failed cases.

Independent tests cover literal identity, signed quarter-turns, basis/mixed-axis half-turns, and a cyclic 120-degree rotation. All largest-component branches include nonzero mixed components, including Y/Z candidates larger than both W and X. Rodrigues' formula supplies independent matrix and vector answers for arbitrary axes, positive/negative and near-half-turn angles. Tests verify unit orientation modulo sign, q/-q matrix equivalence, round trips, orthogonality, determinant +1, quaternion/matrix vector-rotation agreement, Matrix3/Matrix4 extraction, input preservation, independent result storage, and exact float32 exports. Nonunit and zero toMatrix cases characterize unchanged arithmetic. Only Python 3.12 was executed. Two existing quaternion SyntaxWarnings appeared when the edited module was recompiled; their normalization/power checks remain outside scope.

Remaining failures are **37 confirmed defect cases and six contract-question cases**. Confirmed IDs: C02; N01/N02/N03; Q03–Q06/Q08/Q09; R01/R02; E01–E08. Unresolved test contracts cover V01 dimension/empty semantics, V04 ownership, Q10 accuracy, Q11 return semantics, and R03 homogeneous promotion. Detailed identities and counts remain in the result file.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/quaternion.py` | Largest-component selection, branch comparisons, and ctypes refresh |
| `tests/test_quaternion.py` | Convert Q02/Q07 xfails |
| `tests/test_quaternion_matrix.py` | Independent rotation, branch, round-trip, and export regressions |
| `audit/CONVENTIONS.md` | Conversion domain, signs, ownership, and exports |
| `audit/PHASE2-DECISIONS.md` | Record conversion scope without expanding QD06/QD07 policies |
| `audit/COMPATIBILITY.md` | Numerical and ctypes compatibility implications |
| `audit/PHASE2E1.md` | Verification and changed-file summary |
| `audit/phase2e1-test-results.json` | Machine-readable results and remaining failures |
