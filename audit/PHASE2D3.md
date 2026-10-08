# Projection and unprojection correctness

Branch `fix/phase2d3-projection` starts from master `ec0c4e2934caca08e5ea0eff824b588a68d3a909`, the merge of PR #15. This change addresses P01/P02 and the projection input/depth questions in QD03.

## Changes

Projection and unprojection both use row-vector composition modelview*projection. Project divides clip coordinates by W, maps X/Y to the viewport and NDC depth to [0,1], and returns exactly three components. Unproject reverses this mapping and divides the inverse-transformed coordinates by W. Both functions accept Matrix4 wrappers and raw 4×4 nested lists, including mixed representations, without changing their public signatures.

Project requires explicit Vector4 input. The viewport has lower-left origin and upward Y; depth is not clamped. Zero-W and singular-inverse behaviors follow [CONVENTIONS.md](CONVENTIONS.md#projection-and-unprojection). [COMPATIBILITY.md](COMPATIBILITY.md#projection-and-unprojection) describes corrected numerical behavior and migration from compensating operand order or depth mapping. Production changes are confined to these two functions in `gem/matrix.py`.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Unchanged master: **618 passed, 54 strict xfails**. Under `--runxfail`, both P01 input forms and the three-component return case raise TypeError, while P02's noncommuting inverse yields X=-5 instead of 1. All initial 80 new regressions fail before correction.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/projection-results.xml
```

Full suite: **753 cases, 703 passed, 0 failed, 50 strict xfails**, exit 0, with no XPASS or skips. Converts two defect cases (P01/P02), resolves two P01 input/depth question cases, and adds 81 passing cases.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/projection-runxfail.xml
```

This run produced **703 passed, 50 failed**, expected exit 1. JUnit case identities verify that only the four relevant cases leave the expected-failure set and that every unrelated case still fails when exposed. `git diff --check` passes. [Machine-readable results](phase2d3-test-results.json) list all remaining cases and counts by finding and module.

Independent tests use literal windows for identity, asymmetric orthographic bounds, vertical/horizontal perspective, near/far depth, offset nonsquare viewports, and a quarter-turn-plus-translation modelview that does not commute with projection. Sixteen seeded camera cases use analytic eye-space mappings as separate project and unproject oracles, then check round trips. Six general projective matrices use an independent Fraction Gauss-Jordan inverse. Other cases exercise homogeneous scaling, explicit w=0, no Vector3 promotion, input/ctypes preservation, out-of-range depth, zero W, fresh sentinel storage, and singular matrices that project successfully but cannot unproject. Only Python 3.12 was executed.

Remaining failures are **44 confirmed defect cases and six contract-question cases**. Confirmed IDs: C02; M05/M06; N01/N02/N03; Q02–Q09; R01/R02; E01–E08. Unresolved test contracts concern V01 dimension/empty semantics, V04 ownership, Q10 accuracy, Q11 return semantics, and R03 homogeneous promotion. Their detailed identities and counts remain in the result file.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/matrix.py` | Projection, inverse composition, viewport/depth mapping, and matrix input forms |
| `tests/test_matrix.py` | Convert P02 and resolved P01 input/depth cases |
| `tests/test_wiki_contracts.py` | Convert P01 three-component result regression |
| `tests/test_projection.py` | Independent projection, inverse, round-trip, and error regressions |
| `audit/CONVENTIONS.md` | Projection/unprojection contract and transformation cross-reference |
| `audit/PHASE2-DECISIONS.md` | Record QD03's settled ordinary-input contract |
| `audit/COMPATIBILITY.md` | Input, depth, composition, and error implications |
| `audit/PHASE2D3.md` | Verification and changed-file summary |
| `audit/phase2d3-test-results.json` | Machine-readable results and remaining failures |
