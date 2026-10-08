# Plane correctness

Branch `fix/phase2d1-planes` starts from master `d799470c08434ddda061f2fa9cfdb9d06b83ef81`, after PR #13 and the subsequent LICENSE update. This change addresses G01–G04. Historical wiki and source review distinguish definite geometry defects from the representation and offset questions recorded in QD05.

## Changes

Coefficient construction stores scalar `a,b,c,d` and a matching raw normal. Three-point construction stores a unit normal and `d=-normal.dot(point)`. Normalization divides all four coefficients by the original normal length and synchronizes `.normal`. Polygon normal calculation wraps the final vertex to the first.

`bestFitD` retains signed D for `n·p=D`; polygon construction uses coefficient `d=-D`. Public signatures and return behavior are preserved. The field-type correction in `fromPoints`, normal-scale convention, and corrected normalized offset can affect existing callers; migration details are in [COMPATIBILITY.md](COMPATIBILITY.md#plane-construction-and-normalization). [CONVENTIONS.md](CONVENTIONS.md#plane-representation) documents representation, winding, polygon construction, and unresolved error-policy limits.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. The unchanged master baseline produced **440 passed, 61 strict xfails**. All five original G01–G04 cases failed under `--runxfail`; the initial 29 new geometry cases also failed before implementation.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest -q -p no:cacheprovider \
  -o junit_family=legacy --junitxml=/tmp/planes-results.xml
```

Full suite: **545 cases, 489 passed, 56 strict xfails**, exit 0, with no failures, XPASS, or skips. Converted five confirmed xfails (G01: one, G02: two, G03: one, G04: one) and added 44 passing cases. Remaining xfails are 48 other defect cases and eight contract questions.

```bash
/workspace/.venvs/pyGameMath/bin/python -m pytest --runxfail -q --tb=no \
  -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/planes-runxfail.xml
```

This run produced **489 passed, 56 failed**, expected exit 1. JUnit case comparison verifies that every unrelated xfail remains and matches the exposed failures exactly. `git diff --check` passes. [Machine-readable results](phase2d1-test-results.json) record counts by module and finding.

Regression coverage includes positive and negative offsets; raw coefficient normals; all three construction points; positive/negative coefficient scaling; normalization, signed distance, and side classification; clone/flip consistency; polygon closure, cyclic shifts, and reversed winding; signed bestFitD weighting; input independence; and existing degenerate exceptions. Twenty seeded planes use the analytic equation `z=A*x+B*y+C` as an independent normal, incidence, and distance oracle. Four polygon known answers include an oblique plane. Only Python 3.12 was executed.

## Changed files

| File | Purpose |
| --- | --- |
| `gem/plane.py` | Coefficient, construction, normalization, and polygon wrapping fixes |
| `tests/test_plane_ray.py` | Convert five G01–G04 xfails to ordinary regressions |
| `tests/test_planes.py` | Forty-four independent geometry and compatibility cases |
| `audit/CONVENTIONS.md` | Plane representation and polygon construction conventions |
| `audit/PHASE2-DECISIONS.md` | Resolve representation portion of QD05 |
| `audit/COMPATIBILITY.md` | Field, scale, offset, and migration implications |
| `audit/PHASE2D1.md` | Change summary and verification |
| `audit/phase2d1-test-results.json` | Machine-readable suite results |
