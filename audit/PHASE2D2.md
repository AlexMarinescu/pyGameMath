# Matrix translation and vector transformations

Branch `fix/phase2d2-transformations` starts from master `6fcb144a9af4047318d2daf6c86f15e9b37298d6`, the merge of PR #14. This change addresses M04 and V03.

## Changes

In-place Matrix4 translation with Vector3 now constructs a Matrix4 helper, preserving postmultiplication and ctypes synchronization. Vector transformation computes the row-vector product without adding the final row twice. Local implicit w=1 promotion supports affine positions; explicit homogeneous inputs compute output w normally. No perspective divide occurs. Returning and in-place mutation behavior and raw list arguments are preserved.

[CONVENTIONS.md](CONVENTIONS.md#vector-transformations) records the transformation portion of QD09. [COMPATIBILITY.md](COMPATIBILITY.md#matrix-translation-and-vector-transformations) describes changed numerical behavior and migration from compensating transposes or offsets. Projection/unprojection, pivot rotation, shear, and broader shape/error policies remain separate findings.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Unchanged master: **489 passed, 56 strict xfails**. Both original M04/V03 tests fail under `--runxfail`: Vector3 in-place translation raises IndexError, and identity transformation produces `[2,3,5]` instead of `[2,3,4]`. The first 99 new regressions produced **85 failed, 14 passed** before implementation.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/transforms-results.xml
```

Full suite: **672 cases, 618 passed, 0 failed, 54 strict xfails**, exit 0, with no XPASS or skips. Two converted xfails and 127 new cases pass. One existing SyntaxWarning concerns zero-normalization's `is not 0` check (N01), outside this change.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/transforms-runxfail.xml
```

This run produced **618 passed, 54 failed**, expected exit 1. JUnit case identities verify that the remaining xfails equal the baseline minus exactly M04/V03, and every remaining case fails when exposed. `git diff --check` passes. [Machine-readable results](phase2d2-test-results.json) record all remaining failed cases and counts by finding and module.

Independent regression coverage uses literal matrix and angle answers, scaled row basis images, exact integer linearity, analytic affine mappings, and deterministic seeded inputs. Cases cover identity; signed translation on identity, affine, and nonsymmetric nonaffine receivers; returning/in-place ownership and ctypes float32 snapshots; dimensions 1, 2, 3, 4, and 8; known rotations; noncommuting compositions; implicit positions; w=0 directions and arbitrary w; and projective matrices with zero, nonunit, and position-dependent output w. General multiplication's Vector3/Matrix4 failure and existing Matrix2/Matrix3 translation behavior are preserved. Only Python 3.12 was executed.

Remaining failures are **46 confirmed defect cases and eight contract questions**. Confirmed IDs: C02; M05/M06; N01/N02/N03; P01/P02; Q02–Q09; R01/R02; E01–E08. Contract cases cover V01 dimension/empty semantics, V04 ownership, P01 input/depth policy, Q10 accuracy, Q11 return semantics, and R03 homogeneous promotion. Their detailed identities and counts are in the result file.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/matrix.py` | Correct M04 helper dimensions |
| `gem/vector.py` | Row-vector transformation and homogeneous rules |
| `tests/test_matrix.py` | Convert M04 xfail |
| `tests/test_vector_common.py` | Convert V03 xfail |
| `tests/test_transformations.py` | Independent transformation and compatibility regressions |
| `audit/CONVENTIONS.md` | Transformation conventions and dimensional boundaries |
| `audit/PHASE2-DECISIONS.md` | Resolve transformation portion of QD09 |
| `audit/COMPATIBILITY.md` | Numerical behavior and migration implications |
| `audit/PHASE2D2.md` | Verification and changed-file summary |
| `audit/phase2d2-test-results.json` | Machine-readable results and remaining failures |
