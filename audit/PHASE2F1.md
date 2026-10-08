# Ray copying and rigid transforms

R01 duplicates exact stored geometry/distance with independent Vector
storage. R02 rotates start/direction using the Hamilton conjugate sandwich,
keeping Vector fields. R03 applies pure Matrix4 translations with local
w=1 positions and w=0 directions, preserving direction and distance.

Base: master `0d192d92d6058cc6c819c1af5338fe2d535a2ff6` (PR #22).
Branch: `fix/phase2f1-rays`.

## Historical evidence and compatibility

Source `4253839:gem/ray.py` has the same defective copying/transform code;
the archived ray wiki contains only “Coming soon.” Existing matrix rotation
applies to start and direction, establishing a coordinate-origin pivot.
No production consumer establishes intersection state semantics.

Constructor ownership/in-place direction normalization and Matrix3 rotation
remain unchanged, verified by AST comparison. Duplication bypasses
construction so normalized direction does not overwrite stored distance or
state. The misspelled `roateUsingMatrix` and None transform returns remain.
No general Matrix*Vector promotion or quaternion algorithm change is made.

`.end` remains a zero-initialized intersection placeholder/state, copied
exactly by duplication and unchanged by transforms. It does not represent
`start+dir*distance`. There is no way to distinguish a placeholder from a
valid hit at zero or infer validity from a nonzero value. Hit validity,
ownership and transformation policy remain a separate API design question;
no intersection semantics are introduced. Tests calculate geometric tips
independently without using `.end` as a geometric oracle.

See [public API and examples](../docs/RAYS.md),
[conventions](CONVENTIONS.md#ray-copying-and-rigid-transforms) and
[compatibility](COMPATIBILITY.md#phase-2f-1-ray-geometry).

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Baseline: 1,282 passed,
33 expected failures. Both R01 cases, R02 and the R03 promotion question
produce four failures with --runxfail. Before correction, the initial 31 added cases
produce 25 failures and six passes. Four additional Matrix3 intersection-state
checks complete coverage of zero/nonzero state for every transform.

Independent cases cover exact duplication of non-normalized/manual state,
caller ownership, identity/quarter/half/arbitrary-axis rotations, signed
rotations, literal Matrix3 mappings, positive/negative/zero translation,
noncommuting compositions, direction/distance preservation, zero/nonzero
intersection-state storage, independent geometric tips, input preservation,
ctypes snapshot preservation and unchanged general multiplication dimensions.

Full suite: **1,321 passed, 0 failed, 29 expected failures** (1,350 cases).
Only two R01 cases, one R02 case and the R03 promotion-question marker are
converted. The exact unrelated set is preserved. Diagnostic --runxfail:
1,321 passed and exactly those 29 failures (26 defects, three contracts).
[Detailed results](phase2f1-test-results.json) record all remaining identities.

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/ray-final.xml
python -m pytest --runxfail -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/ray-runxfail.xml
```

## Changed files

| File | Change |
| --- | --- |
| `gem/ray.py` | Exact isolated duplication, quaternion rotation, local translation |
| `tests/test_plane_ray.py` | Convert R01/R02/R03 regressions |
| `tests/test_api_edges.py` | Convert R01 aliasing regression |
| `tests/test_rays.py` | 35 independent geometry/ownership regressions |
| `docs/RAYS.md` | Public ray semantics, examples and unresolved hit state |
| `audit/CONVENTIONS.md` | Copy/rigid-transform contract |
| `audit/COMPATIBILITY.md` | Numerical and ownership implications |
| `audit/PHASE2-DECISIONS.md` | QD02/QD04 boundaries and end ambiguity |
| `audit/PHASE2F1.md` | Historical evidence and verification |
| `audit/phase2f1-test-results.json` | Counts and remaining case identities |
