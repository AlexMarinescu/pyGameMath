# Phase 2C angle conversion and vector refraction

Branch: `fix/phase2c-angles-refraction`, based on latest master `1f657b17755bc4238979eb77e2abfcaaf90c4b01`, the merge of PR #12. Prior audit reports, compatibility, mathematical conventions, decision checklist, and historical wiki evidence were reviewed before implementation. Scope stops at C01 and V02.

## Implementation and approved contract

The angle helpers now multiply by full-precision math.pi conversion factors, preserving original names, signatures and legacy keyword parameter names. No other API's angle unit or algorithm changes.

Refraction uses the standard vector equation with IOR squared, and local Vector*scalar operands instead of unsupported scalar*Vector. The historical wiki confirms the argument order but does not settle ratio direction or unit/orientation requirements. The ambiguity was reported before implementation, and the user explicitly approved IOR=n1/n2, normalized incident and normal vectors, the normal opposing incidence and pointing into the incident medium, and the historical zero Vector result for total internal reflection. No automatic normalization or normal flipping is introduced.

See [CONVENTIONS.md](CONVENTIONS.md#phase-2c-current-angle-and-refraction-conventions) for the equation, air/glass examples, interface direction, and total-internal-reflection sentinel, and [COMPATIBILITY.md](COMPATIBILITY.md#phase-2c-c01-v02) for changed answers and migration implications. The refraction portion of QD10 is recorded as approved; its viewport question and all other unresolved policies remain deferred. Critical-angle rounding uses the existing `k<0` branch without a new tolerance/clamp.

## Validation

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Unchanged merged baseline: **356 passed / 65 strict xfails**, exit 0. Before corrections, focused `--runxfail` execution of all four C01/V02 cases reproduced **4 failed**, exit 1.

Full suite after corrections:

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/phase2c-results.xml
```

**501 cases: 440 passed, 61 strict xfails, no failures, XPASS, or skips; exit 0.** Four original xfails are now ordinary passes (C01: two, V02: two), and 80 new cases pass. Remaining xfails comprise **53 other confirmed-defect cases and eight contract questions**. Baseline/final JUnit case comparison verifies that exactly C01/V02 cases leave the expected-failure set; every unrelated marked case remains.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy \
  --junitxml=/tmp/phase2c-runxfail.xml
```

**440 passed / 61 failed**, expected exit 1. Failures match the remaining expected-failure set exactly. `git diff --check` passes. [Machine-readable summary](phase2c-test-results.json) records counts by module and finding. Only Python 3.12 was executed; interpreter support metadata is unchanged.

New test coverage:

- Thirteen known positive/negative/zero angles in both directions, nine round-trip values including small and large ordinary finite inputs, 20 seeded standard-library comparison/odd-symmetry cases, and legacy keyword compatibility.
- Five normal-incidence ratios; four oblique known answers including a non-axis-aligned normal; explicit air-to-glass and glass-to-air at normal and 30-degree incidence; total internal reflection beyond the critical angle at three ratios including glass-to-air; and a boundary input whose discriminant evaluates to zero.
- Twenty independently constructed 3D normal/tangent geometries derive transmitted direction from Snell's law and trigonometry. Check direction, unit length, tangential Snell relation, orientation, reciprocal reverse travel, and input preservation. No library normalization, rotation, or conversion helper supplies the geometry oracle.

The original total-internal-reflection regression remains an ordinary pass. Existing mixed-unit rotation/projection and all prior matrix/quaternion regressions remain unchanged. No NumPy, Cython, or compiled dependency is introduced.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/common.py` | C01 conversion formulas and signature-preserving docstrings |
| `gem/vector.py` | V02 squared ratio, local operand order, and contract docstring |
| `tests/test_vector_common.py` | Convert four C01/V02 xfails to ordinary regressions |
| `tests/test_angles_refraction.py` | Eighty known-answer/property/ownership cases |
| `audit/CONVENTIONS.md` | Current angle/refraction contracts and medium-transition examples |
| `audit/PHASE2-DECISIONS.md` | Record explicit approval of refraction portion of QD10 |
| `audit/COMPATIBILITY.md` | Changed results, retained signatures, and migration guidance |
| `audit/PHASE2C.md` | Scope, approval provenance, tests, and complete file summary |
| `audit/phase2c-test-results.json` | Machine-readable validation summary |

Historical reports, snapshots, and raw baseline results are preserved. No unrelated fixes, further phase work, merge into master, or PyPI publication are included.
