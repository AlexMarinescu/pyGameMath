# Pivot rotation and shear correctness

Branch `fix/phase2d4-pivot-shear` starts from master `efa68429f0a314a2cdd3289643cbcaf9a9636722`, the merge of PR #16. This change addresses M05/M06 and records QD09's pivot and shear mappings.

## Changes

`rotate2` computes the row-vector offset `pivot-pivot*R`, keeping nonzero pivots stationary. XY3 and XY4 write shear factors into the Z column; XY3 no longer indexes past its rows, and XY4 preserves homogeneous W. The existing YZ/XZ mappings and class rotation dispatch are preserved. Signatures, row-major storage, matrix product order, and ctypes refresh behavior are unchanged.

Historical source `4253839` contains the same defective formulas as the pre-fix master. The wiki names rotation and shear operations without defining pivot semantics or displacement coordinates. [CONVENTIONS.md](CONVENTIONS.md#pivot-rotation-and-shear) defines those mappings. [COMPATIBILITY.md](COMPATIBILITY.md#pivot-rotation-and-shear) documents the changed XY4 behavior and migration from pivot-offset workarounds. Production mathematics changes only in rotate2, shearXY3, and shearXY4.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Unchanged master: **703 passed, 50 strict xfails**. Both original cases fail under `--runxfail`: XY3 raises IndexError, and a 90-degree pivot transform maps `(2,3)` to `(4,6)`. The 95 new regressions produced **47 failed, 48 passed** before correction.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/pivot-results.xml
```

Full suite: **848 cases, 800 passed, 0 failed, 48 strict xfails**, exit 0, with no XPASS or skips. Two converted M05/M06 cases and 95 new cases pass.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/pivot-runxfail.xml
```

This run produced **800 passed, 48 failed**, expected exit 1. JUnit case identities verify that the remaining xfails equal the baseline minus exactly M05/M06, and every unrelated case still fails when exposed. `git diff --check` passes. [Machine-readable results](phase2d4-test-results.json) list all remaining cases and counts by finding and module.

Independent regression coverage includes literal quarter/half/eighth-turn answers; positive/negative angles and pivots; stationary pivots and directions; analytic relative-coordinate rotation, distance preservation, and inverse composition; unchanged Matrix2 and axis-angle dispatch; all shear planes with positive/negative/zero factors; basis-image and coordinate oracles on identity/nonsymmetric receivers; arbitrary W preservation; returning/in-place methods, caller storage, and exact float32 exports; noncommuting shear/shear and shear/translation compositions; and determinant/inverse properties with deterministic inputs. Only Python 3.12 was executed.

Remaining failures are **42 confirmed defect cases and six contract-question cases**. Confirmed IDs: C02; N01/N02/N03; Q02–Q09; R01/R02; E01–E08. Test contracts concern V01 dimension/empty semantics, V04 ownership, Q10 accuracy, Q11 return semantics, and R03 homogeneous promotion. Detailed identities and counts remain in the result file.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/matrix.py` | Pivot offsets, XY3/XY4 shear entries, and mapping docstrings |
| `tests/test_matrix.py` | Convert M05/M06 xfails |
| `tests/test_pivot_shear.py` | Independent pivot, shear, composition, and compatibility regressions |
| `audit/CONVENTIONS.md` | Pivot/shear mappings, dimensions, and rotation dispatch |
| `audit/PHASE2-DECISIONS.md` | Record QD09 pivot/shear decisions |
| `audit/COMPATIBILITY.md` | Changed pivot and XY4 numerical behavior |
| `audit/PHASE2D4.md` | Verification and changed-file summary |
| `audit/phase2d4-test-results.json` | Machine-readable results and remaining failures |
