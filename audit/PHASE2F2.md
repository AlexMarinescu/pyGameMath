# Finite numerical robustness

N01 restores direct zero-vector and identity-quaternion normalization.
N02 uses scaled hypot norms and normalization. N03 adds exact power-of-two
scaling around the existing 3x3/4x4 cofactor kernels, preserving public
signatures, row-vector conventions and storage.

Base: master `2cbd899a47546b8c40f8244f57799969c86b0d87` (PR #23).
Branch: `fix/phase2f2-numerical-robustness`.

## Evidence and implementation

Historical `4253839` initializes zero/identity fallback results but its
`length is not 0` guard incorrectly enters division. Phase 1B identifies
source fallback evidence; the wiki defines no broader exception policy.
Core callers previously relied on this exception for invalid axes/directions.
Narrow exact-zero guards preserve those errors without changing ray geometry,
plane equations, rotation mathematics, interpolation or experimental code.

Repeated unscaled hypot rounds the norm of three minimum subnormals to one
minimum subnormal, rather than two. Scaling before hypot avoids that loss;
normalization independently scales components to preserve direction even
when a finite norm would overflow. Legacy nonfinite-input arithmetic remains.

Arbitrary row-max division introduces rounding that can turn singular
[[1,2,3],[4,5,6],[7,8,9]] into a nonzero floating determinant. Binary power
scaling avoids that additional rounding. Floating cofactor arithmetic can
still miss singularity (including duplicated rows at extreme scales), so
an exact integer determinant check clears binary denominators from each
row of the represented binary64 input. This is a zero check, not a tolerance
or a different inverse algorithm. Existing cofactor kernels remain unchanged,
verified by AST comparison after renaming. 2x2 inverse and public determinant
functions are also unchanged. Severe conditioning and extreme relative
scales remain limitations of cofactor inversion; no accuracy promise is
made for unrepresentable results. Signed infinity handles rescaling overflow.

See [conventions](CONVENTIONS.md#numerical-robustness) and
[compatibility](COMPATIBILITY.md#phase-2f-2-numerical-robustness).

## Accuracy and performance

Known nonsymmetric triangular answers at uniform scales 1e-300 through
1e300 differ by at most 1.66e-16 relative error in the benchmark.
Independent regression tolerances are 2e-15 for norms/unit directions and
3e-14 for nonsymmetric inverse/identity products. Subnormal rounding and
float32 export limits are tested explicitly.

Final ordinary-scale benchmark, median process-CPU microseconds over five
runs of 3,000 calls:

| Size | Baseline raw | Final raw | Baseline wrapper | Final wrapper |
| --- | --- | --- | --- | --- |
| 3 | 1.367 | 9.355 | 3.380 | 12.011 |
| 4 | 5.975 | 21.851 | 8.129 | 23.832 |

The exact singularity check increases cost beyond the earlier scaling-only
prototype. Raw kernels cost about 6.8x/3.7x; wrappers about 3.6x/2.9x.
These environment-specific measurements describe the correctness trade-off.
[Benchmark data](phase2f2-inversion-benchmarks.json) are reproducible with:

```sh
python audit/benchmark-numerical-inversion.py
```

## Testing

CPython 3.12.14, pytest 9.1.1, six 1.17.0. Baseline: 1,321 passed and
29 expected failures. All seven original N01/N02/N03 cases reproduce their
failures. Running the 55 new cases against an isolated baseline archive
produces 37 failures and 18 passes.

Decimal (800-digit) norms/directions and exact Fraction inverses provide
independent references. Tests cover zeros and ordinary/unit inputs, extremes,
mixed scales, minimum subnormals, overflow, nonsymmetric/diagonal matrices,
both inverse multiplication orders, singular matrices, degenerate callers,
returning/in-place ownership, legacy nonfinite norms and ctypes exports.

Full suite: **1,383 passed, 0 failed, 22 expected failures** (1,405 cases).
Only N01 (two), N02 (three) and N03 (two) markers are converted. The exact
unrelated set is unchanged. Diagnostic --runxfail: 1,383 passes and exactly
those 22 failures (19 defects, three contracts). [Detailed results](phase2f2-test-results.json)
list every remaining identity. N03's abs=0-only comparison is made explicitly
relative (1e-14), since binary64 cofactor results do not promise bit equality.

```sh
python -m pytest -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/numeric-final.xml
python -m pytest --runxfail -q -p no:cacheprovider -o junit_family=legacy --junitxml=/tmp/numeric-runxfail.xml
```

## Changed files

| Files | Change |
| --- | --- |
| `gem/vector.py`, `gem/quaternion.py` | Stable normalization/norms; zero fallbacks and axis guards |
| `gem/matrix.py` | Scaled inverse wrappers; exact singular checks; axis/lookAt guards |
| `gem/plane.py`, `gem/ray.py`, `gem/common.py` | Narrow exact-zero guards |
| `tests/test_vector_common.py`, `tests/test_quaternion.py`, `tests/test_matrix.py` | Convert only numerical xfails |
| `tests/test_numerical_robustness.py` | 55 independent cases |
| `docs/QUATERNIONS.md` | Current direct normalization contract |
| `audit/CONVENTIONS.md`, `audit/COMPATIBILITY.md`, `audit/PHASE2-DECISIONS.md` | Numerical boundaries and compatibility |
| `audit/PHASE2F2.md`, `audit/phase2f2-test-results.json` | Verification and remaining failures |
| `audit/benchmark-numerical-inversion.py`, `audit/phase2f2-inversion-benchmarks.json` | Reproducible final benchmark |
