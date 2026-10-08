# Quaternion axes and arithmetic

Branch `fix/phase2e2-quaternion-operations` starts from master `c018af877e2e44cc9efefd144667d9ccf3fc2f29`, the merge of PR #18. This change addresses Q06/Q08/Q09 and the axis-helper ownership portion of QD02.

## Changes

Axis-angle helpers normalize temporary Vector/list axes and preserve caller storage. The arbitrary-axis helper retains its existing normalized-axis Quaternion sandwich result, leaving Q11 unresolved. In-place quaternion/vector multiplication uses Vector component storage. Python 3 division aliases the retained legacy float-only methods. Quaternion order, Hamilton multiplication, public signatures, and established degree/radian units are unchanged.

Historical source `4253839` shows the same broken list branches and Vector-axis mutation; the quaternion wiki remains a placeholder. Baseline characterization confirms that both helpers replace caller Vector storage with unit components. [CONVENTIONS.md](CONVENTIONS.md#quaternion-axes-and-arithmetic) records temporary-axis ownership and arithmetic protocols. [COMPATIBILITY.md](COMPATIBILITY.md#quaternion-axes-and-arithmetic) describes migration from incidental normalization. Production changes are confined to these operations in `gem/quaternion.py`.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Unchanged master: **872 passed, 43 strict xfails**. All four original regressions fail under `--runxfail`: two list-axis AttributeErrors, a missing Vector `.data` AttributeError, and Python 3 division TypeError. All 118 new cases fail before correction.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/quat-ops-results.xml
```

Full suite: **1033 cases, 994 passed, 0 failed, 39 strict xfails**, exit 0, with no XPASS or skips. Converts two Q06 cases, one Q08 case, and one Q09 case; adds 118 passing regressions.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/quat-ops-runxfail.xml
```

This run produced **994 passed, 39 failed**, expected exit 1. JUnit case comparison verifies that only Q06/Q08/Q09 leave the expected-failure set; Q11 and every unrelated case remain and fail when exposed. `git diff --check` passes. [Machine-readable results](phase2e2-test-results.json) record module/finding counts and all remaining failed cases.

Independent regression coverage includes signed standard and arbitrary nonunit axes, list/Vector inputs, caller storage identity, zero/unsupported axes, unchanged legacy arbitrary-axis output, degree/radian differences, positive axis rescaling, Hamilton scalar-dot/vector-cross known answers, deterministic products, receiver identity, legacy/Python 3 division, positive/negative/zero operands, unsupported types, float subclasses, and failure atomicity. Only Python 3.12 was executed. Two existing normalization/power SyntaxWarnings appeared on recompilation and remain outside scope.

Remaining failures are **33 confirmed defect cases and six contract-question cases**. Confirmed IDs: C02; N01/N02/N03; Q03/Q04/Q05; R01/R02; E01–E08. Test contracts concern V01 dimension/empty semantics, V04 ownership, Q10 accuracy, Q11 return semantics, and R03 homogeneous promotion. Detailed identities and counts remain in the result file.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/quaternion.py` | Temporary axis handling, Vector product storage, and Python 3 division |
| `tests/test_quaternion.py` | Convert Q06/Q08/Q09 xfails |
| `tests/test_quaternion_operations.py` | Independent axis, arithmetic, ownership, and error regressions |
| `audit/CONVENTIONS.md` | Axis and arithmetic protocols |
| `audit/PHASE2-DECISIONS.md` | Record axis-helper ownership boundary |
| `audit/COMPATIBILITY.md` | Ownership, product, and division implications |
| `audit/PHASE2E2.md` | Verification and changed-file summary |
| `audit/phase2e2-test-results.json` | Machine-readable results and remaining failures |
