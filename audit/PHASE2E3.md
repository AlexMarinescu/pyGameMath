# Quaternion powers and logarithms

Q03 powers now scale the input imaginary axis using the principal angle,
including fractional and negative powers. Q04 logarithms return zero scalar
and independent storage. The supported domain is unit quaternions, with
explicit identity, negative-identity and zero handling; see
[conventions](CONVENTIONS.md#quaternion-powers-and-logarithms) and
[compatibility](COMPATIBILITY.md#phase-2e-3-quaternion-powers-and-logarithms).

Base: master `649b6d1808175291b8badacfd3453f7c15a4681f` (PR #19).
Branch: `fix/phase2e3-quaternion-powers`.

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Baseline: 994 passed,
39 expected failures. Running the three original Q03/Q04 regressions with
`--runxfail` produced three failures. Before correction, the first 130 new
cases produced 94 failures and 36 passes; one overflow regression was added
after implementation.

The 131 added cases cover independently repeated Hamilton products, known
axis logarithms, integer/fractional/negative powers, pure-real identities,
zero rejection, sign-sensitive branches, tiny imaginary components down to
subnormal scale, fresh result storage, and finite extreme exponent output.

Full suite: **1,128 passed, 0 failed, 36 expected failures** (1,164 cases).
The remaining set exactly matches baseline minus the two Q03 and one Q04
cases. A diagnostic `--runxfail` run produced exactly those same 36 failures
and 1,128 passes: 30 confirmed-defect cases and six unresolved-contract
cases. Detailed case identities are in [test results](phase2e3-test-results.json).

Commands (using `/workspace/.venvs/pyGameMath/bin/python`):

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/powers-final.xml
python -m pytest --runxfail -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/powers-runxfail.xml
```

## Changed files

| File | Change |
| --- | --- |
| `gem/quaternion.py` | Unit power/log formulas and singular-axis handling |
| `tests/test_quaternion.py` | Remove only Q03/Q04 strict defect markers |
| `tests/test_quaternion_powers.py` | Independent numerical and ownership regressions |
| `audit/CONVENTIONS.md` | Principal-angle and unit-domain conventions |
| `audit/COMPATIBILITY.md` | Numerical, storage and error implications |
| `audit/PHASE2-DECISIONS.md` | Power/log scope within QD06/QD07 |
| `audit/PHASE2E3.md` | Verification and change summary |
| `audit/phase2e3-test-results.json` | Counts and exact remaining failures |
